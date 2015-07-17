import os, sys, logging, re
import requests
import json
import argparse
import datetime
import subprocess
import xml.dom.minidom

token = None
headers = None
lane_tags = None
flowcell_lane_cache = dict()
flowcell_contenttype = None
base_api_url = None

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger('upload_data.py')

script_options = {
    "base_api_url": None,
    "basedir": os.getcwd(),
    "quiet": False,
    "debug": False,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-a", "--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.")

    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required.")

    return parser

class ClusterMonitor(object):

    QSTAT_COMMAND = ["qstat", "-xml"]

    def __init__(self, api_url, token):
        self.api_url = api_url
        self.token = token
        self.headers = {'Authorization': "Token %s" % token}

    def run(self):
        currently_running = self.parse_jobnames()
        marked_running = self.lims_currently_processing()

        # What alignments are currently running but not marked yet?
        mark_alignments = currently_running - marked_running
        if mark_alignments:
            log.info("Marking alignments as processing: %s" % str(mark_alignments))
            [self.update_processing_status(align_id, True) for align_id in mark_alignments]

        # What alignments are currently marked but not running?
        finished_alignments = marked_running - currently_running
        if finished_alignments:
            log.info("Alignments no longer processing: %s" % str(finished_alignments))
            [self.update_processing_status(align_id, False) for align_id in finished_alignments]

    def update_processing_status(self, align_id, processing=True):

        patch_url = "%s/flowcell_lane_alignment/%d/" % (self.api_url, align_id)

        update_result = requests.patch(patch_url, headers = self.headers, data={'currently_processing': processing})

        if not update_result.ok:
            log.critical("Could not update alignment %d: %s" % (align_id, str(update_result)))
            return False

        return True

    def lims_currently_processing(self):

        fetch_url = "%s/flowcell_lane_alignment/?currently_processing=True" % self.api_url

        fetch_results = requests.get(fetch_url, headers = self.headers)
        lims_process_align_ids = set()

        if fetch_results.ok:
            results = fetch_results.json()
            for result in results['results']:
                lims_process_align_ids.add(result['id'])
            log.info("Currently marked as processing on LIMS: %s" % str(lims_process_align_ids))
            return lims_process_align_ids
        else:
            log.critical("Could not execute api query: %s" % fetch_url)
            sys.exit(1)

    def parse_jobnames(self):
        qstat_xml = subprocess.check_output(self.QSTAT_COMMAND, stderr=subprocess.STDOUT)
        dom=xml.dom.minidom.parseString(qstat_xml)

        jobnames = dom.getElementsByTagName('JB_name')

        alignments = set()

        for job in jobnames:
            jobname = job.childNodes[0].data
            log.debug(jobname)
            match = re.search("ALIGN#(\d+)", jobname)
            if match:
                alignments.add(int(match.group(1)))

        log.info("Alignment IDs: %s" % alignments)
        return alignments

def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the default logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)
        # Make this a little less noisy by default
        requests_log = logging.getLogger("requests.packages.urllib3.connectionpool")
        requests_log.setLevel(logging.WARN)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
        log.debug("Using LIMS API endpoint: %s from environment" % api_url)
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
        log.debug("Using LIMS API endpoint: %s from options" % api_url)
    else:
        sys.stderr.write("Could not find LIMS API URL.\n")
        sys.exit(1)


    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        sys.stderr.write("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    monitor = ClusterMonitor(api_url, token)

    monitor.run()

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
