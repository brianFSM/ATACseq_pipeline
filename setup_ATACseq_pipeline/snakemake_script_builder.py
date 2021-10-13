'''
Creates Snakefiles for ATACseq pipeline with default settings
'''
# from os import mkdir
import argparse
import pathlib
# import time
# import string
from jinja2 import Environment, FileSystemLoader


def init_argparse():
    '''
    function docstring goes here
    '''
    parser = argparse.ArgumentParser(
        prog="python auto_run_mutect2.py",
        usage="%(prog)s --sample_base 'sample_base'",
        description="Creates, and submits to slurm, a script for \
            running mutect2 based on the number of tumor/normal \
            samples"
        )
    parser.add_argument(
        "-v", "--version", action="version",
        version="{parser.prog} version 1.0.0"
        )
    parser.add_argument('-s','--sample_list', nargs='+', help='<Required> Set flag', required=True)
    # Use like:
    # python arg.py -l 1234 2345 3456 4567
    parser.add_argument(
        "-r", "--reference", action='store',
        default="hg38",
        help="The reference genome."
        )
    parser.add_argument(
        "-m", "--mode", action='store',
        default="PE",
        help="Sequence mode, single-end (SE) or paired-end (PE)"
    )
    return parser


def render_output(project_path,
                  template_path,
                  template_filename,
                  output_filename,
                  argument_dict):

    '''
    function docstring goes here
    '''
    # Environment for jinja
    env = Environment(
        loader=FileSystemLoader(template_path))

    # Set up some paths
    template_path = str(template_path) + "/" + template_filename
    output_path = str(project_path) + "/" + output_filename

    # Open and read our template
    _template = env.get_template(template_filename)
    template_str = open(template_path, mode='r').read()


    # First grab the tumor bam filenames and turn them into a mutect2
    # argument
    end_line_string = "\n"
    samples_string = ""

    for sample in argument_dict["samples"]:
        print(f"sample = {sample}")
        fq_path="data/reads/" + sample + "_R1.fq.gz"
        if argument_dict["mode"] == "PE":
            fq_path += ", data/reads/" + sample + "_R2.fq.gz"

        samples_string += "      " + sample + ": " + \
                    "[" + fq_path + "]" + end_line_string
        print(f"samples_string: {samples_string}")

    # Process the template into our new output
    new_output_str = _template.render(reference=argument_dict["reference"],
                                      samples=samples_string)

    with open(output_path, "w") as new_file:
        new_file.write(new_output_str)



def main():
    '''
    docstring goes here
    '''
    # parse the arguments and store them
    parser = init_argparse()
    args = parser.parse_args()

    reference = args.reference
    samples = args.sample_list[0].split()
    mode = args.mode

    # Setup the directories for input and output
    current_wd = pathlib.Path().absolute()
    # parent_path = str(current_wd.parent)
    # new_dir_path = parent_path + "/" + project_title + "_" + \
        # str(time.time()) + "/"

    arg_dict = {"samples": samples,  "reference": reference, "mode": mode}
    template_path = \
    "/home/bwp9287/pipelines/ATACseq_pipeline/testDir/setup_ATACseq_pipeline/templates"
    output_file = "configs/config_samples.yaml"

    render_output(current_wd, template_path, "config_samples_template.yaml",
                  output_file, arg_dict)


if __name__ == "__main__":
    main()
