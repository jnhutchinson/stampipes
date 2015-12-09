import json
import os
import sys
import argparse
import logging
import requests
import subprocess

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')

script_options = {
    "quiet": False,
    "debug": False,
    "base_api_url": None,
    "token": None,
    "aggregation_ids": [],
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "overwrite": False,
    "script_name": "run.bash",
    "qsub_prefix": ".agg",
    "dry_run": False,
    "aggregation_base_directory": None,
    "aggregation_directory": None,
    "script_template": None,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-a", "--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.  Required if not in the environment under LIMS_API_URL")
    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required if not in the environment under LIMS_API_TOKEN")
    parser.add_argument("--aggregation_base_directory", dest="aggregation_base_directory",
        help="The base directory to put aggregations in.  Can get from environment AGGREGATIONS variable")

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="Append commands to run this alignment to this file.")
    parser.add_argument("--overwrite", dest="overwrite", action="store_true",
        help="Create a new outfile instead of appending commands.")
    parser.add_argument("--script_template", dest="script_template",
        help="The script template to use.")
    parser.add_argument("--aggregation_directory", dest="aggregation_directory",
        help="The directory for the aggregation.  Will deduce if not given.")
    parser.add_argument("-b", "--script_basename", dest="script_name",
        help="Name of the script that goes after the sample name.")

    parser.add_argument("--tag", dest="tag",
        help="Run for alignments tagged here.")    
    parser.add_argument("--aggregation", dest="aggregation_ids", type=int, action="append",
        help="Run for these aggregations (can be used more than once).")

    parser.add_argument("--qsub-prefix", dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
        help="Take no action, only print messages.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


class ProcessSetUp(object):

    def __init__(self, args, api_url, token, aggregation_base_directory): 

        self.token = token
        self.api_url = api_url
        self.qsub_scriptname = args.script_name
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.dry_run = args.dry_run
        self.aggregation_base_directory = aggregation_base_directory
        self.aggregation_directory = args.aggregation_directory
        self.headers = {'Authorization': "Token %s" % self.token}
        self.script_template = args.script_template
        self.overwrite = args.overwrite

    def api_single_result(self, url_addition=None, url=None):

        if url_addition:
           url = "%s/%s" % (self.api_url, url_addition)

        request = requests.get(url, headers=self.headers)

        if request.ok:
            logging.debug(request.json())
            return request.json()
        else:
            logging.error("Could not get data from %s" % url)
            logging.error(request)
            return None

    def api_list_result(self, url_addition=None, url=None):

        more = True
        results = []

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        while more:

            logging.debug("Fetching more results for query %s" % url)

            request = requests.get(url, headers=self.headers)

            if not request.ok:
                logging.error(request)
                return None
            more_results = request.json()
            results.extend(more_results["results"])
            if more_results["next"]:
                url = more_results["next"]
            else:
                more = False

        return results

    def get_aggregation_info(self, aggregation_id):

        results = self.api_single_result("aggregation/%d" % aggregation_id)

        if not results:
            logging.error("Could not find information for aggregation %d" % aggregation_id)
            return None

        return results

    def set_aggregation_folder(self, aggregation_info, library_info):

        url = "directory/?purpose__slug=all-alignments-bam&filetype__slug=bam&content_type=%d&object_id=%d" % (aggregation_info["object_content_type"], aggregation_info["id"])
        results = self.api_list_result(url)

        if len(results) > 1:
            logging.error("Found %d folders for aggregation %d, require 1" % (len(results), aggregation_info["id"]))

        if len(results) == 0:

            if not self.aggregation_base_directory:
                logging.critical("Connot proceed, no directory set and no base aggregation directory given.")
                sys.exit(1)

            purpose_url = "file_purpose/?slug=all-alignments-bam"
            results = self.api_list_result(purpose_url)

            if not results:
                logging.critical("Cannot find file purpose all-alignment-bam")
                sys.exit(1)

            file_purpose = results[0]["url"]

            path = os.path.join(self.aggregation_base_directory, "LN%d" % library_info["number"], "aggregation-%d" % aggregation_info["id"])

            logging.info("Setting aggregation folder to %s" % path)

            data = {
                "content_type": "%s/content_type/%d/" % (self.api_url, aggregation_info["object_content_type"]),
                "object_id": aggregation_info["id"],
                "path": path,
                "purpose": file_purpose,
            }

            new_result = requests.post("%s/directory" % self.api_url, headers = self.headers, data = data)

            if not new_result.ok:
                logging.critical(new_result)
                logging.critical("Could not upload new aggregation folder path to LIMS: %s" % json.dumps(data))
                sys.exit(1)

            return path

        return results[0]["path"]

    def get_aggregation_lanes(self, aggregation_id):
        results = self.api_list_result("aggregation_lane/?aggregation=%d&include=True" % aggregation_id)

        if not results:
            logging.error("Could not find lanes for aggregation %d" % aggregation_id)
            return []

        return results

    def get_lane_alignments_file(self, aggregation_id, alignment_id):

        results = self.api_list_result("file/?purpose__slug=all-alignments-bam&filetype__slug=bam&content_type=47&object_id=%d" % alignment_id)

        if len(results) != 1:
            logging.error("Found %d files for alignment %d, require 1 (Aggregation %d)" % (len(results), alignment_id, aggregation_id))
            logging.error(results)
            return None
        file_info = results[0]

        return (file_info["path"], file_info["md5sum"])

    def get_library_info(self, aggregation_info):
        library_info = self.api_single_result(url=aggregation_info["library"])
        if not library_info:
            logging.critical("Cannot proceed without library!  Could not get info from %s (Aggregation %d)" % (aggregation_info["library"], aggregation_info["id"]))
            sys.exit(1)
        return library_info

    def get_genome_index(self, aggregation_info):
        genome_info = self.api_single_result(url=aggregation_info["genome_index"])
        if not genome_info:
            logging.critical("Could not get genome info! (Aggregation %d)" % aggregation_info["id"])
            sys.exit(1)
        return genome_info

    def get_genome_index_location(self, aggregation_id, aggregation_lanes):
        included = None
        for aggregation_lane in aggregation_lanes:
            if aggregation_lane["include"]:
               included = aggregation_lane
               break

        alignment = self.api_single_result(url=aggregation_lane["alignment"])

        if not alignment:
            logging.critical("Was not able to fetch alignment %s! (Aggregation %d)" % (aggregation_lane["alignment"], aggregation_id))
            sys.exit(1)

        genome_location = self.api_single_result(url=alignment["genome_index_location"])
        if not genome_location:
            logging.critical("Could not get genome location from alignment %d! (Aggregation %d)" % (included["id"], aggregation_id))
            sys.exit(1)

        return os.path.join(genome_location["base_path"], genome_location["directory"], genome_location["filename"])

    def get_script_template(self, aggregation_id, process_template_url, script_template=None):

        if script_template:
            logging.info("Using script template %s" % script_template)
            return open(script_template, 'r').read()

        if not process_template_url:
            logging.critical("No process template for aggregation %d\n" % aggregation_id)
            return None

        logging.info("Getting process template %s" % process_template_url)

        process_template = self.api_single_result(url=process_template_url)

        if not process_template:
            logging.critical("Could not find processing template for %s\n" % process_template_url)
            return None

        script_path = os.path.expandvars(process_template["process_version"]["script_location"])
        return open(script_path, 'r').read()

    def get_example_flowcell(self, aggregation_id, aggregation_lanes):
        included = None
        for aggregation_lane in aggregation_lanes:
            if aggregation_lane["include"]:
                included = aggregation_lane
                break

        lane = self.api_single_result(url=aggregation_lane["lane"])

        if not lane:
            logging.critical("Was not able to fetch lane %s (Aggregation %d)" % (aggregation_lane["lane"], aggregation_id))
            sys.exit(1)

        flowcell = self.api_single_result(url=lane["flowcell"])
        if not flowcell:
            logging.critical("Could not get flowcell at %d (Aggregation %d)" % (lane["flowcell"], aggregation_id))
            sys.exit(1)

        return flowcell

    def add_script(self, aggregation_id, aggregation_folder, library_number):

        if self.overwrite:
            mode = "w"
        else:
            mode = "a"

        with open(self.outfile, mode) as runfile:
            runfile.write("cd %s && " % aggregation_folder)
            runfile.write("qsub -V -cwd -S /bin/bash -N \"%sLN%d_AGG#%d\" %s\n\n" % (self.qsub_prefix, library_number, aggregation_id, self.qsub_scriptname))

    def setup_tag(self, tag_slug):

        aggregation_tags = self.api_list_result("tagged_object?content_type=126&tag__slug=%s" % tag_slug)

        [self.setup_aggregation(aggregation_tag["object_id"]) for aggregation_tag in aggregation_tags]

    def setup_aggregation(self, aggregation_id):

        aggregation = self.get_aggregation_info(aggregation_id)

        if not aggregation:
            return False

        aggregation_lanes = self.get_aggregation_lanes(aggregation_id)

        if not aggregation_lanes:
            return False

        library_info = self.get_library_info(aggregation)
        aggregation_folder = self.set_aggregation_folder(aggregation, library_info)
        genome_index = self.get_genome_index(aggregation)
        genome_index_location = self.get_genome_index_location(aggregation_id, aggregation_lanes)
        flowcell = self.get_example_flowcell(aggregation_id, aggregation_lanes)

        logging.info("Aggregation %d folder: %s" % (aggregation_id, aggregation_folder))
        logging.debug(aggregation)

        missing = False

        files = []
        for aggregation_lane in aggregation_lanes:
            if not aggregation_lane["include"]:
                logging.info("Not including lane %s (Aggregation %d)" % (aggregation_lane["lane"], aggregation_id))
                continue
            alignment_endpoint = aggregation_lane["alignment"]

            if not alignment_endpoint:
                logging.info("Not including lane %s because no alignment set (Aggregation %d)" % (aggregation_lane["lane"], aggregation_id))

            alignment_id = int(alignment_endpoint.strip("/").split("/")[-1])

            bamfile = self.get_lane_alignments_file(aggregation_id, alignment_id)

            if not bamfile:
                logging.critical("No BAM alignment file for alignment %s for lane %s, skipping (Aggregation %d)" % (alignment_endpoint, aggregation_lane["lane"], aggregation_id))
                missing = True
                continue
            else:
                logging.info(bamfile)
                files.append(bamfile)

        if missing: return False

        script_contents = self.get_script_template(aggregation_id, aggregation["aggregation_process_template"], self.script_template)

        if not script_contents:
            logging.critical("No script contents")
            return

        os.makedirs(aggregation_folder, exist_ok=True)

        file_record = open("%s/bamfiles.txt" % aggregation_folder, "w")
        file_record.write("\n".join(["\t".join(bamfile) for bamfile in files]))
        file_record.close()

        script_file = os.path.join(aggregation_folder, self.qsub_scriptname)
        logging.info("Creating script file %s" % script_file)

        script = open(script_file, "w")
        script.write("export AGGREGATION_ID=%d\n" % aggregation_id)
        script.write("export LIBRARY=%d\n" % library_info["number"])
        script.write("export LIBRARY_NAME=LN%d\n" % library_info["number"])
        script.write("export BAM_FILES=\"%s\"\n" % " ".join([bamfile[0] for bamfile in files]))
        script.write("export GENOME=%s\n" % genome_index["label"])
        script.write("export GENOME_INDEX=%s\n" % genome_index_location)
        script.write("export AGGREGATION_FOLDER=%s\n" % aggregation_folder)
        script.write("export READ_LENGTH=%d\n" % flowcell["read_length"])
        if aggregation["umi"]:
            script.write("export UMI=True\n")
        else:
            script.write("unset UMI\n")
        if flowcell["paired_end"]:
            script.write("export PAIRED=True\n")
        else:
            script.write("unset PAIRED\n")
        script.write("\n")

        script.write(script_contents)

        script.close()

        self.add_script(aggregation_id, aggregation_folder, library_info["number"])

def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
        logging.getLogger("requests").setLevel(logging.WARNING)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)
        logging.getLogger("requests").setLevel(logging.WARNING)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
    else:
        logging.error("Could not find LIMS API URL.\n")
        sys.exit(1)

    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        logging.error("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    if not poptions.aggregation_base_directory and "AGGREGATIONS" in os.environ:
        aggregation_base_dir = os.environ["AGGREGATIONS"]
    elif poptions.aggregation_base_directory:
        aggregation_base_dir = poptions.aggregation_base_directory
    else:
        aggregation_base_dir = None

    process = ProcessSetUp(poptions, api_url, token, aggregation_base_dir)

    for aggregation_id in poptions.aggregation_ids:
        process.setup_aggregation(aggregation_id)

    if poptions.tag:
        process.setup_tag(poptions.tag)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
