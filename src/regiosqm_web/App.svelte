<script>
    import ModalBase from './components/ModalBase.svelte'
    import Modal from './components/Modal.svelte'
    import {openModal, closeModal, closeAllModals, modals} from './stores/modals.js'

    // import Editor from './components/Editor.svelte'
    import Navigation from './components/Navigation.svelte'
    export let appTitle = 'Regioselectitiy App'

    const handleClick = () => {
        openModal(Modal, {title: 'Alert', message: 'This is an alert'})
    }
    const handleClick2 = () => {
        openModal(Modal, {title: 'Alert', message: 'This another alert'}, true)
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
                <p><a href="https://fu">Read more...</a></p>
                <button on:click={handleClick}>Open Modal</button>
                <button on:click={handleClick2}>Other Modal</button>
            </div>
            <div class="lg:w-2/3 md:w-1/2 bg-gray-300 rounded-lg sm:mr-10 relative">
                <!-- <Editor /> -->
            </div>
        </div>
    </section>
</div>

{#if $modals.length > 0}
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
</style>
