from argparse import ArgumentParser
from pathlib import Path
from subprocess import Popen, PIPE
from os import getenv
from datastructures import model_path
import logging
import inlist_templates
import json

logging.basicConfig(level='INFO')

pre_parser = ArgumentParser(add_help=False)

pre_parser.add_argument('--template', type=str)

args, _ = pre_parser.parse_known_args()

parser = ArgumentParser()

parser.add_argument('--template', type=str, required=True)

parser.add_argument('--model-dir',
                    type=Path,
                    default=model_path)

parser.add_argument('--mesa-object', type=str, required=True)
parser.add_argument('--mesa-index', type=str, required=True)
parser.add_argument('--name', type=str, required=False)

if args.template:
    template = inlist_templates.get_template('Gyre', args.template)
    template.apply_args(parser)

args = parser.parse_args()

if not args.model_dir.is_dir():
    logging.error("Model directory %s does not exist", args.model_dir)
    exit(1)

if not args.model_dir.exists():
    logging.error("Model directory %s is not a directory", args.model_dir)
    exit(1)

if not (gyre_dir := getenv("GYRE_DIR")):
    logging.error("Missing GYRE_DIR environment variable (see https://gyre.readthedocs.io/en/stable/user-guide/quick-start.html )")
    exit(1)

mesa_model = (args.model_dir.joinpath(args.mesa_object)
                            .joinpath('LOGS')
                            .joinpath(f'profile{args.mesa_index}.data.GYRE'))
model = args.model_dir.joinpath(args.mesa_object, 'gyre', args.mesa_index, template.get_suggested_model_name(args) + (("_" + args.name) if args.name else ''))
model.mkdir(parents=True)
inlist = model.joinpath('inlist')

with inlist.open('w') as f:
    f.write(template.render(model, mesa_model, args))

gyre = Popen([f'{getenv("GYRE_DIR")}/bin/gyre', inlist], text=True, stdout=PIPE)

gyre_version = gyre.stdout.readline().split()[1]

json.dump({
    'type': 'Gyre',
    'version': gyre_version,
    'template': {
        'name': template.name,
        'version': template.version
    },
    'base': {
        'model': args.mesa_object,
        'index': args.mesa_index
    },
    'properties': template.get_properties(args)
}, model.joinpath('properties.json').open('w'))

for line in gyre.stdout:
    logging.info("%s", line[:-1])
