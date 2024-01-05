<script>
    // import ModalBase from './components/ModalBase.svelte'

    // TODO refer to global stores and components? global name space?

    import Modal from '$components/Modal.svelte'
    import {openModal, closeModal, modals} from '$stores/modals.js'
    // import {rdkit, jquery, chemdoodle} from './stores/libs.js'

    import ChemdoodleEditor from './features/ChemdoodleEditor/index.svelte'

    import {Grid, Row, Column} from 'carbon-components-svelte'
    import {Button, ButtonSet} from 'carbon-components-svelte'
</script>

<svelte:head>
    <title>RegioSQM</title>
    <link rel="stylesheet" href="/fontawesome/css/all.css" />
    <link rel="stylesheet" href="/fonts/font.css" />
    <meta name="google" content="notranslate" />
    <meta name="mobile-web-app-capable" content="yes" />
    <meta name="theme-color" content="#ffffff" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no" />
</svelte:head>

<Grid fullWidth narrow>
    <Row>
        <Column>
            <h1 class="logo">RegioSQM</h1>
        </Column>
        <Column>
            <ButtonSet>
                <Button kind="ghost">Usage</Button>
                <Button kind="ghost">FAQ</Button>
                <Button kind="ghost">GitHub</Button>
            </ButtonSet>
        </Column>
    </Row>
    <Row>
        <Column>
            <!-- Entry text -->
            <h2 class="center bold">Where would regioselective attack happen?</h2>
            <span class="h-12 block" />
        </Column>
    </Row>
    <Row>
        <Column>
            <ChemdoodleEditor />
        </Column>
    </Row>
    <Row>
        <Column>
            <span class="h-12 block" />
        </Column>
    </Row>
    <Row>
        <Column>
            <span class="text-center h-12 block leading-12"> Have feedback? Tweet at me bro or open issue </span>
        </Column>
    </Row>
</Grid>

<style global lang="postcss">
    @import 'tailwindcss/base';
    @import 'tailwindcss/components';
    @import 'tailwindcss/utilities';
    @import 'carbon-components-svelte/css/g10.css';
    /* @import 'carbon-components-svelte/css/white.css'; */

    /* help guides */
    .bx--grid,
    .bx--row,
    .bx--col {
        box-shadow: inset 0px 0px 0px 1px rgba(0, 0, 0, 0.5);
    }

    body {
        font-family: Nunito, serif;
        /* background: #f9f8f2; */
        background: #f1f1f1;
        color:#343433;
    }

    .bx--btn {
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

    .logo {
        font-size: 25px;
        line-height: 48px;
        font-weight: 700;
        letter-spacing: 0;
    }

    .leading-12 {
        line-height: 3rem;
    }

    p {
        @apply my-5;
    }
</style>
