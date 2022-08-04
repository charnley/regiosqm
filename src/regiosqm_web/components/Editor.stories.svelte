<script>
    import {Meta, Template, Story} from '@storybook/addon-svelte-csf'
    import Editor, {chemdoodleGetMol, chemdoodleSetMol} from './Editor.svelte'

    import {rdkit, jquery, chemdoodle} from '../stores/libs.js'

    const onRdkitLoaded = () => {
        // if RDKit is not undefined
        initRDKitModule().then((instance) => {
            rdkit.update(() => instance)
        })
    }

    const onChemdoodleLoaded = () => {
        isChemdoodleCoreLoaded = true
    }

    const onChemdoodleLoaded2 = () => {
        chemdoodle.update(() => ChemDoodle)
        jquery.update(() => ChemDoodle.lib.jQuery)
    }

    const checkAllGlobals = () => {
        if ($rdkit == null) return false
        if ($jquery == null) return false
        if ($chemdoodle == null) return false

        pageReady = true
    }

    let isChemdoodleCoreLoaded = false
    let pageReady = false

    $: $rdkit, checkAllGlobals()
    $: $jquery, checkAllGlobals()
    $: $chemdoodle, checkAllGlobals()
</script>

<svelte:head>
    <script src="/rdkit/RDKit_minimal.js" on:load={onRdkitLoaded}></script>
    <script src="/chemdoodleweb/ChemDoodleWeb-unpacked.js" on:load={onChemdoodleLoaded}></script>
    {#if isChemdoodleCoreLoaded}
        <script src="/chemdoodleweb/ChemDoodleWeb-uis-unpacked.js" on:load={onChemdoodleLoaded2}></script>
    {/if}
</svelte:head>
<Meta title="2D Editor" component={Editor} />

<Template let:args>
    <Editor />
</Template>

<Story name="Primary" />
