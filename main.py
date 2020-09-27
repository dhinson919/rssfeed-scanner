#!/usr/bin/env python

"""
Checks configured RSS feeds for disease-specific channel item activity.

Usage:  $ ./main.py [-d | --debug (debug logging)]

Sample output:
    $ ./main.py
    INFO:root:Begin checking disease RSS feed activity; inactivity period = 14 days
    INFO:root:[Wolman Disease] has been inactive
    WARNING:root:[Cron's Disease] is not a recognized MeSH disease descriptor name; skipping
    WARNING:root:[Crohn's Disease] is not a recognized MeSH disease descriptor name; skipping
    INFO:root:End checking disease RSS feed activity

Comments:

The stated requirement is
    "Given a tuple of Disease and RSS feeds, determine which disease had no activity for a given number of days."

To that end I've made the following assumptions:

1.  A tuple will have the form (disease-name, RSS-feed-1 [, RSS-feed-n]...)
2.  The workflow will involve multiple such tuples in a given cycle.
3.  The use of the term "Disease" is significant, as opposed to "a topic such as Disease".  Therefore the expectation
    is that the workflow will only query RSS feeds concerning diseases.  Since it is unclear whether the supplied
    disease-names will be guaranteed to represent diseases I decided to confirm it myself.  I do this by querying
    the NIH's Medical Subject Headings database to ensure that the disease-name matches scientifically canonical
    disease names.  If a provided disease-name does not match a lookup against that database then I do not attempt to
    query RSS feeds about it.

TODO:

I consider this to be a prototype for how such a workflow might be implemented at a high level.  For a production-ready
system I believe it would need (at least) the following improvements that are outside the time scope of this work:

1.  Automated unit tests.  Some manual testing has been performed but it was not rigorous.
2.  Parallelization.  The workflow performs only serialized access to RSS feeds.  A concurrent worker model
    implemented either internally or externalized to more appropriate infrastructure would be very desirable.
3.  Robustness.  There is basic error detection and reporting but no mechanisms for retrying failed communication
    attempts beyond what underlying modules automatically provide.
4.  Scheduling.  There's no provision for automating the scheduling of the workflow.
5.  Logging.  Only rudimentary logging has been implemented.
6.  Metrics.  There are no performance or accounting metrics being collected.

"""

__author__ = "David Hinson (27-SEP-2020)"
__email__ = "dhinson@sjmail.us"
__status__ = "Prototype"

import collections
import getopt
import logging
import re
import sys
import xml.etree.ElementTree as ET
from functools import reduce

import requests

# Template for NIH MeSH service tree query:
Treenum_query_template = '''
PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
SELECT ?treeNum
FROM <http://id.nlm.nih.gov/mesh>
WHERE {{ mesh:{descriptor_id} meshv:treeNumber ?treeNum }}
ORDER BY ?treeNum
'''

# Regex for parsing MeSH base tree IDs:
Base_tree_regex = re.compile('http://id.nlm.nih.gov/mesh/([^.]+)')
# Regex for paring MeSH descriptor IDs:
Descriptor_id_regex = re.compile('http://id.nlm.nih.gov/mesh/(\D\d+)')

# Templates for RSS feed service API query parameters indexed by feed URL.  As new RSS feed services are identified and
# added appropriate templates must be added here.
FeedQueryParams = {
    'https://clinicaltrials.gov/ct2/results/rss.xml': 'lup_d={inactivity_period_days}&cond={disease}&count={limit}'
}


def fetch_mesh_descriptor(descriptor_label):
    """
    Fetches the Medical Subject Heading (MeSH) descriptor for the provided descriptor label from the NIH MeSH service.
    :param descriptor_label: Speculative name of a descriptor label
    :return: MeSH descriptor object if the descriptor label is know, otherwise empty dict
    """
    result = {}

    api_query_params = {'label': descriptor_label, 'match': 'exact', 'limit': 10}
    response = None
    try:
        response = requests.get('https://id.nlm.nih.gov/mesh/lookup/descriptor', params=api_query_params, timeout=10)
    except Exception as exc:
        logging.error(exc)
    if response and response.status_code == requests.codes.ok:
        logging.debug('Call to https://id.nlm.nih.gov/mesh/lookup/descriptor took {seconds} seconds'
                      .format(seconds=response.elapsed))

        # If a JSON object is found then return it directly:
        response_json = response.json()
        if response_json:
            result = response_json[0]
    return result


def is_mesh_disease(mesh_descriptor):
    """
    Queries the NIH MeSH RDF service to test if the provided MeSH descriptor represents a known disease.
    :param mesh_descriptor: MeSH Descriptor of a speculative disease
    :return: True if the descriptor represents a known disease
    """
    result = False

    match = Descriptor_id_regex.search(mesh_descriptor['resource'])
    if match:
        descriptor_id = match.group(1)
        sparql_query = Treenum_query_template.format(descriptor_id=descriptor_id)
        api_query_params = {'query': sparql_query}
        response = None
        try:
            response = requests.get('https://id.nlm.nih.gov/mesh/sparql', params=api_query_params, timeout=10)
        except Exception as exc:
            logging.error(exc)
        if response and response.status_code == requests.codes.ok:
            logging.debug('Call to https://id.nlm.nih.gov/mesh/sparql took {seconds} seconds'
                          .format(seconds=response.elapsed))

            # Find the first descriptor tree number in the XML response and test if it is in the MeSH Diseases tree:
            content_stringed = response.content.decode()
            xml_root = ET.fromstring(content_stringed)
            treenum_uri_nodes = xml_root.findall(
                ".//"
                "{http://www.w3.org/2005/sparql-results#}binding[@name='treeNum']/"
                "{http://www.w3.org/2005/sparql-results#}uri")
            if bool(treenum_uri_nodes):
                treenum_id = treenum_uri_nodes[0].text
                match = Base_tree_regex.search(treenum_id)
                if match:
                    base_tree = match.group(1)
                    if base_tree.startswith('C'):
                        result = True
    return result


def count_channel_items(feed_url, disease_name, inactivity_period_days):
    """
    Counts the disease-related channel items of the specified RSS feed. 
    :param feed_url: URL of the RSS feed
    :param disease_name: Name of the disease
    :param inactivity_period_days: Inactivity period threshold in days
    :return: A tuple containing the qualifying channel item count and the value certainty. If the count is zero and
    certainty is True then there were no channel items found. If the count is zero and the certainty is False then
    the actual count is unknown due to an unsuccessful query.
    """
    count = 0
    certainty = False
    if feed_url not in FeedQueryParams:
        logging.error('No query params found for feed [{feed_url}] while checking [{disease_name}]'
                      .format(feed_url=feed_url, disease_name=disease_name))
    else:
        api_query_params = FeedQueryParams[feed_url].format(
            disease=disease_name, inactivity_period_days=inactivity_period_days, limit=1000)
        response = None
        try:
            response = requests.get(feed_url, params=api_query_params, timeout=10)
        except Exception as exc:
            logging.error(exc)
        if response and response.status_code == requests.codes.ok:
            logging.debug('Call to {feed_url} took {seconds} seconds'
                          .format(feed_url=feed_url, seconds=response.elapsed))

            # Search XML response for existence of channel items:
            content_stringed = response.content.decode()
            xml_root = ET.fromstring(content_stringed)
            item_pubdate_nodes = xml_root.findall('./channel/item')
            count = len(item_pubdate_nodes)
            certainty = True
    return count, certainty


# A custom iterator capable of being preempted early.
class SmartIterator(collections.Iterable):
    class __Inner:
        __underlying_iterator = None
        __stop = False

        def __init__(self, iterable):
            self.__underlying_iterator = iterable.__iter__()

        def __next__(self):
            if self.__stop:
                raise StopIteration
            else:
                return self.__underlying_iterator.__next__()

        def stop(self):
            self.__stop = True

    __iterable = None
    _inner = None

    def __init__(self, iterable):
        self.__iterable = iterable

    def __iter__(self):
        self._inner = SmartIterator.__Inner(self.__iterable)
        return self._inner


# SmartIterator for use by feed iteration.
class FeedsIterator(SmartIterator):
    def __call__(self, channel_items):
        count, _ = channel_items
        if count > 0:
            self._inner.stop()
        return channel_items


def check_disease_feeds(disease_name, feed_urls, inactivity_period_days):
    """
    Check provided RSS feeds for disease-related channel items.
    :param disease_name: Name of disease topic
    :param feed_urls: List of feed URLs
    :param inactivity_period_days: Inactivity period threshold in days
    """
    # Map RSS feed activity to list of (channel item count, certainty) tuples.  We augment the list comprehension
    # with a "smart" iterator so that evaluation of subsequent feeds can be avoided should a preceding feed be found
    # to contain activity, satisfying the per-disease test.
    feeder = FeedsIterator(feed_urls)
    activity = [feeder(count_channel_items(feed_url, disease_name, inactivity_period_days)) for feed_url in feeder]
    # Sum list of (channel item count, certainty) tuples. Note that uncertainty is propagated:
    channel_items_total, is_certain = reduce((lambda acc, arg: (acc[0] + arg[0], acc[1] and arg[1])), activity)
    # Log inactivity results:
    if channel_items_total == 0:
        if is_certain:
            logging.info("[{disease_name}] has been inactive".format(disease_name=disease_name))
        else:
            logging.info("Inactivity for [{disease_name}] could not be determined with certainty"
                         .format(disease_name=disease_name))


def check_disease_activity(work_unit, inactivity_period_days):
    """
    Evaluate a disease activity check unit of work.
    :param work_unit: Disease activity check unit of work
    :param inactivity_period_days: Inactivity period threshold in days
    """
    disease_name, feed_urls = work_unit
    # Attempt to match speculative disease name with MeSH descriptor of a recognized disease:
    mesh_descriptor = fetch_mesh_descriptor(disease_name)
    if not bool(mesh_descriptor) or not is_mesh_disease(mesh_descriptor):
        logging.warning("[{disease_name}] is not a recognized MeSH disease descriptor name; skipping"
                        .format(disease_name=disease_name))
    else:
        # Proceed with checking feeds assuming MeSH disease name as working name:
        mesh_disease_name = mesh_descriptor['label']
        check_disease_feeds(mesh_disease_name, feed_urls, inactivity_period_days)


def log_disease_inactivity(work_units, inactivity_period_days):
    """
    Main workflow entry point.
    :param work_units: List of work unit tuples containing speculative disease names and associated RSS feed URLS
    :param inactivity_period_days: Inactivity period threshold in days
    """
    logging.info('Begin checking disease RSS feed activity; inactivity period = {period} days'
                 .format(period=inactivity_period_days))
    for work_unit in work_units:
        check_disease_activity(work_unit, inactivity_period_days)
    logging.info('End checking disease RSS feed activity')


def main(argv):
    log_level = logging.INFO

    try:
        opts, args = getopt.getopt(argv[1:], 'd')
    except getopt.GetoptError as exc:
        sys.exit(exc)

    for opt, arg in opts:
        if opt in ('-d', '--debug'):
            log_level = logging.DEBUG

    logging.basicConfig(level=log_level)

    sample_work_units = [
        # Case with currently active channel items:
        ('Alzheimer Disease', ['https://clinicaltrials.gov/ct2/results/rss.xml']),
        # Case with currently inactive channel items:
        ('Wolman Disease', ['https://clinicaltrials.gov/ct2/results/rss.xml']),
        # Case with incorrect disease name:
        ("Cron's Disease", ['https://clinicaltrials.gov/ct2/results/rss.xml']),
        # Case with common but non-MeSH-preferred disease name:
        ("Crohn's Disease", ['https://clinicaltrials.gov/ct2/results/rss.xml']),
        # Case with currently active channel items involving more than one feed:
        ('Crohn Disease', ['https://clinicaltrials.gov/ct2/results/rss.xml', 'https://clinicaltrials.gov/ct2/results/rss.xml']),
    ]

    log_disease_inactivity(sample_work_units, inactivity_period_days=14)


if __name__ == '__main__':
    main(sys.argv)
