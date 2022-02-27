
# Regiosqm

Lorem ipsum dolor sit amet, officia excepteur ex fugiat reprehenderit enim
labore culpa sint ad nisi Lorem pariatur mollit ex esse exercitation amet. Nisi
anim cupidatat excepteur officia. Reprehenderit nostrud nostrud ipsum Lorem est
aliquip amet voluptate voluptate dolor minim nulla est proident. Nostrud
officia pariatur ut officia. Sit irure elit esse ea nulla sunt ex occaecat
reprehenderit commodo officia dolor Lorem duis laboris cupidatat officia
voluptate. Culpa proident adipisicing id nulla nisi laboris ex in Lorem sunt
duis officia eiusmod. Aliqua reprehenderit commodo ex non excepteur duis sunt
velit enim. Voluptate laboris sint cupidatat ullamco ut ea consectetur et est
culpa et culpa duis.

    - https://doi.org/10.1039/C7SC04156J
    - https://doi.org/10.1186/s13321-021-00490-7
    - https://github.com/jensengroup/regiosqm
    - https://github.com/NicolaiRee/RegioSQM20


## Usage

Can be used as

    - Python Library
    - Python API/computation server
    - Web application

## Dependencies

You need to install conda and Node

## Setup

As per all other repo's I play with, there is a Makefile for quick setup. Just run

    make

and it will setup all the dependencies and run initial build. The following
sections are more in depth.

### Setup Svelite frontend

The JavaScript is a single-page-application (SPA), which means you just need to
compile the JavaScript once and make it avaiable to the user.

    npm install
    npm run build

### Setup Python backend

The backend is based on `ppqm` and `rdkit`.

    conda env create -p ./env

and then you need to setup `ppqm`, which does not have a conda package as this
writing. Clone the repo down and install it with pip (in the proper environment).

    git clone https://github.com/ppqm/ppqm ppqm.git
    cd ppqm.git
    ../env/bin/python -m pip install -e .

## Development

You would want to run `watch` with npm, so it builds

## Deployment

- [ ] TODO guide. For ngix and apache2

## TODO

- [X] Chemdoodle integration / jquery clicks
- [X] PostCSS on generated elements
- [X] AJAX Requests to cactus
- [X] Read result into chemdoodle from cactus
- [x] Integrate rdkit, for smiles convertion and more
- [x] Wait/loader while all external libraries are loading
- [x] popup / modal / userfeedback
- [X] Setup storybook for Svelte components (with global tailwind and local fonts)

- [X] Python3 and ppqm implementation of RegioSQM2020
- [X] Python3 and ppqm implementation of RegioSQM2018
- [ ] Interface to RegioML

- [X] Python Flask server
- [X] Use Flask to serve Svelte Single page application
- [X] Flask SQLAlchemy interface
- [ ] REST API For submitting
- [ ] REST API For checking status
- [ ] REST API For fetching results
- [ ] Slurm-like queuing system interface
- [ ] Only allow API requests from active frontend (SSO Cookie?)

- [ ] About page
- [ ] Result page

- [ ] Read up on web workers for heavy tasks (load the full lib per worker)


## Rdkit Note


JavaScript functions on a Mol obj

    add_hs: ƒ Mol$add_hs()
    condense_abbreviations: ƒ ()
    draw_to_canvas: ƒ Mol$draw_to_canvas(arg0, arg1, arg2)
    draw_to_canvas_with_highlights: ƒ Mol$draw_to_canvas_with_highlights(arg0, arg1)
    draw_to_canvas_with_offset: ƒ Mol$draw_to_canvas_with_offset(arg0, arg1, arg2, arg3, arg4)
    generate_aligned_coords: ƒ ()
    get_aromatic_form: ƒ Mol$get_aromatic_form()
    get_cxsmarts: ƒ Mol$get_cxsmarts()
    get_cxsmiles: ƒ Mol$get_cxsmiles()
    get_descriptors: ƒ Mol$get_descriptors()
    get_inchi: ƒ Mol$get_inchi()
    get_json: ƒ Mol$get_json()
    get_kekule_form: ƒ Mol$get_kekule_form()
    get_molblock: ƒ Mol$get_molblock()
    get_morgan_fp: ƒ ()
    get_morgan_fp_as_uint8array: ƒ ()
    get_new_coords: ƒ ()
    get_smarts: ƒ Mol$get_smarts()
    get_smiles: ƒ Mol$get_smiles()
    get_stereo_tags: ƒ Mol$get_stereo_tags()
    get_substruct_match: ƒ Mol$get_substruct_match(arg0)
    get_substruct_matches: ƒ Mol$get_substruct_matches(arg0)
    get_svg: ƒ ()
    get_svg_with_highlights: ƒ Mol$get_svg_with_highlights(arg0)
    get_v3Kmolblock: ƒ Mol$get_v3Kmolblock()
    is_valid: ƒ Mol$is_valid()
    remove_hs: ƒ Mol$remove_hs()
