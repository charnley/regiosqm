
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


## Setup

### Frontend

    npm install
    npm run dev

### Backend

    conda env create -p ./env
    python -m regiosqm_api --start

## TODO

 - [X] Chemdoodle integration / jquery clicks
 - [X] PostCSS on generated elements
 - [X] AJAX Requests to cactus
 - [X] Read result into chemdoodle from cactus
 - [x] Integrate rdkit, for smiles convertion and more
 - [x] Wait/loader while all external libraries are loading
 - [x] popup / modal / userfeedback

 - [ ] Python backend server
 - [X] Use Flask to serve Svelte Single page application
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
