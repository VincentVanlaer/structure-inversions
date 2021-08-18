from argparse import ArgumentParser
from pathlib import Path
from subprocess import Popen, PIPE
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
parser.add_argument('--work-dir',
                    type=Path)

parser.add_argument('--name', type=str)

if args.template:
    template = inlist_templates.get_template('MESA', args.template)
    template.apply_args(parser)

args = parser.parse_args()

if not args.model_dir.is_dir():
    logging.error("Model directory %s does not exist", args.model_dir)
    exit(1)

if not args.model_dir.exists():
    logging.error("Model directory %s is not a directory", args.model_dir)
    exit(1)

if not args.work_dir:
    args.work_dir = args.model_dir.joinpath('workdir/MESA')

if not args.work_dir.exists():
    logging.error("Work directory %s does not exist", args.model_dir)
    exit(1)

if not args.work_dir.is_dir():
    logging.error("Work directory %s is not a directory", args.model_dir)
    exit(1)

if not args.work_dir.joinpath('star').exists():
    logging.error("Work directory does not look like a MESA work directory, 'star' executable not found")
    exit(1)

name = template.get_suggested_model_name(args)
name = (args.name + "_" + name) if args.name else name
model = args.model_dir.joinpath(name)
model.mkdir()
model.joinpath('LOGS').mkdir()
model.joinpath('photos').mkdir()
inlist = model.joinpath('inlist')


with inlist.open('w') as f:
    f.write(template.render(model, args))

star = Popen([args.work_dir.joinpath('star'), inlist], text=True, stdout=PIPE, cwd=args.work_dir)

mesa_version = star.stdout.readline().split()[1]

json.dump({
    'type': 'MESA',
    'version': mesa_version,
    'template': {
        'name': template.name,
        'version': template.version
    },
    'properties': template.get_properties(args)
}, model.joinpath('properties.json').open('w'))

for line in star.stdout:
    logging.info("%s", line[:-1])
