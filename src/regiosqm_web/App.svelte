<script>
    import ModalBase from './components/ModalBase.svelte'
    import Modal from './components/Modal.svelte'
    import {openModal, closeModal, closeAllModals, modals} from './stores/modals.js'

    import Editor, {chemdoodleGetMol, chemdoodleSetMol} from './components/Editor.svelte'
    import Navigation from './components/Navigation.svelte'
    export let appTitle = 'Regioselectitiy App'

    /**
     *
     * @param query
     * @param resultFormat
     * @param successFunction
     * @param failedFunction
     */
    const requestCactus = (query, resultFormat, successFunction, failedFunction) => {
        // for example https://cactus.nci.nih.gov/chemical/structure/butanol/smiles

        let search

        search = query
        search = search.replace('[', '%5B')
        search = search.replace(']', '%5D')
        search = search.replace('@', '%40')
        search = search.replace('=', '%3D')
        search = search.replace('#', '%23')

        const url = 'https://cactus.nci.nih.gov/chemical/structure/' + search + '/' + resultFormat

        // Tell user to wait
        openModal(Modal, {title: 'Requesting from cactus', message: '', addCancelButton: false})

        fetch(url)
            .then((response) => {
                closeModal()
                if (response.status == 200) {
                    successFunction(response)
                    // response.text().then((text) => {
                    //     let molecule = text
                    //     openModal(Modal, {message: 'Molecule found: ' + molecule}, true)
                    // })
                } else {
                    failedFunction(response)
                    // openModal(Modal, {title: 'Error', message: 'Could not find the molecule on cactus'}, true)
                }
            })
            .catch((error) => {
                failedFunction(undefined)
                console.log(error)
                openModal(
                    Modal,
                    {
                        message: 'Connection problems to cactus.nci.nih.gov. Wait a bit and try again',
                    },
                    true,
                )
                return []
            })
    }

    function smilesToSdf(smi) {
        // var mol = RDKit.Molecule.fromSmiles( smi );
        let mol = RDKit.get_mol(smi)

        // mol.addHs()
        // mol.Embedmolecule3D()
        mol.remove_hs()
        let mol3d = mol.get_molblock()

        return mol3d
    }

    function sdfToSmiles(sdf) {
        // NOTE Seems like removeHs is added automatically
        // NOTE WHY IS RDKIT JAVASCRIPT IN SNAKE_CASE?!
        let mol = RDKit.get_mol(sdf)
        mol.add_hs()
        var smi = mol.get_smiles()
        return smi
    }

    const searchExample = () => {
        let name = 'butane'
        let format = 'smiles'

        let molSmiles = requestCactus(
            name,
            format,
            () => {},
            () => {},
        )

        let sdf = smilesToSdf(molSmiles)
    }

    // searchExample()

    // Await all modules are loaded
    let pageReady = false
    let RDKit

    function onLoaded() {
        initRDKitModule().then((instance) => {
            RDKit = instance
            pageReady = true
        })
    }

    // User experience
    let userSearchQuery

    const handleSearchSubmit = () => {
        requestCactus(
            userSearchQuery,
            'smiles',
            (response) => {
                response.text().then((text) => {
                    let smiles = text
                    let sdf = smilesToSdf(smiles)
                    chemdoodleSetMol(sdf)
                })
            },
            () => {},
        )
    }

    function capitalizeFirstLetter(string) {
        // why is this not a standard javascript function
        return string.charAt(0).toUpperCase() + string.slice(1).toLowerCase()
    }

    const handleWhatClick = () => {
        let mol = chemdoodleGetMol()
        console.log(mol)
        let smiles = sdfToSmiles(mol)
        console.log(smiles)
        requestCactus(
            smiles,
            'iupac_name',
            (response) => {
                response.text().then((text) => {
                    console.log(text)
                    text = capitalizeFirstLetter(text)
                    openModal(Modal, {title: text, message: 'Molecule found on cactus'}, true)
                })
            },
            () => {},
        )
    }
</script>

<svelte:head>
    <title>{appTitle}</title>
    <meta name="theme-color" content="#ffffff" />
    <meta name="google" content="notranslate" />
    <meta name="mobile-web-app-capable" content="yes" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no" />
    <link href="https://fonts.googleapis.com/css?family=Roboto|Fira+Sans:300,400,500" rel="stylesheet" />
    <link href="https://fonts.googleapis.com/css?family=Nunito:400,700,800" rel="stylesheet" />
    <link rel="stylesheet" href="/fontawesome/css/all.css" />
    <script src="RDKit_minimal.js" on:load={onLoaded}></script>
</svelte:head>

<div class="min-h-screen">
    <Navigation />

    <section class="relative mt-20">
        <div class="container flex sm:flex-nowrap flex-wrap max-w-7xl mx-auto px-4 sm:px-6">
            <div class="lg:w-1/3 md:w-1/2 bg-white flex flex-col md:ml-auto w-full">
                <h1 class="text-4xl mb-1">Predict Regioselectitiy</h1>
                <p>
                    of electrophilic aromatic substitution reactions. For citation and method details, please see the
                    RegioSQM/ML papers.
                </p>
                <p>
                    Well, at lest it will be. First POC is just interacting between chemdoodle, rdkit and cactus. E.i.
                    find and get a name for a structure.
                </p>
                <div>
                    <button
                        class="h-10 px-4 text-white rounded-lg bg-red-500 hover:bg-red-600"
                        on:click={handleWhatClick}>What is this Pokemon?!</button>
                </div>
            </div>
            <div class="lg:w-2/3 md:w-1/2 bg-gray-300 rounded-lg sm:mr-10 relative">
                <div class="container flex items-center pt-5 px-5">
                    <div class="relative">
                        <form on:submit|preventDefault={handleSearchSubmit}>
                            <div class="absolute top-4 left-3">
                                <i class="fa fa-search text-gray-400 z-20 hover:text-gray-500" />
                            </div>
                            <input
                                type="text"
                                class="h-14 w-100 pl-10 pr-20 rounded-lg z-0 focus:shadow focus:outline-none"
                                placeholder="Search anything..."
                                bind:value={userSearchQuery} />
                            <div class="absolute top-2 right-2">
                                <button
                                    type="submit"
                                    class="h-10 w-20 text-white rounded-lg bg-red-500 hover:bg-red-600">Search</button>
                            </div>
                        </form>
                    </div>
                </div>

                <Editor />
            </div>
        </div>
    </section>
</div>

{#if $modals.length > 0 || !pageReady}
    <slot name="backdrop" />

    <div class="fixed z-10 inset-0 overflow-y-auto" aria-labelledby="modal-title" role="dialog" aria-modal="true">
        <div class="flex items-end justify-center min-h-screen pt-4 px-4 pb-20 text-center sm:block sm:p-0">
            <!--
Background overlay, show/hide based on modal state.

Entering: "ease-out duration-300"
From: "opacity-0"
To: "opacity-100"
Leaving: "ease-in duration-200"
From: "opacity-100"
To: "opacity-0"
-->
            <div class="fixed inset-0 bg-gray-500 bg-opacity-75 transition-opacity" aria-hidden="true" />

            <!-- This element is to trick the browser into centering the modal contents. -->
            <span class="hidden sm:inline-block sm:align-middle sm:h-screen" aria-hidden="true">&#8203;</span>
            <ModalBase>
                <div slot="backdrop" class="backdrop" on:click={closeModal} />
            </ModalBase>
        </div>
    </div>
{/if}

<style global lang="postcss">
    @tailwind base;
    @tailwind components;
    @tailwind utilities;

    body {
        font-family: Nunito, sans-serif;
        font-weight: bold;
    }

    .backdrop {
        position: fixed;
        top: 0;
        bottom: 0;
        right: 0;
        left: 0;
        background: rgba(0, 0, 0, 0.5);
    }

    p {
        @apply my-5;
    }
</style>
