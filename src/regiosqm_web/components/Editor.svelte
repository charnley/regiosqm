<script>
    import IconBtn from './IconBtn.svelte'

    export let sketcherName = 'sketcherSingle'
    export let atoms = ['H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Br', 'I']

    let clientWidth
    let clientHeight

    let jquery
    let sketcher

    const initializeChemdoodle = () => {
        sketcher = new ChemDoodle.SketcherCanvas(sketcherName, 100, 100, {useServices: false, oneMolecule: true})
        chemdoodleResize(sketcher, [500, 500])
        jquery = ChemDoodle.lib.jQuery
        onWindowResize()
    }

    const chemdoodleResize = (canvas, dim) => {
        var width = dim[0]
        var height = dim[1]
        canvas.resize(width, height)
        setTimeout(function () {
            chemdoodleClick('#sketcherSingle_button_scale_plus')
        }, 50)
    }

    const chemdoodleClick = (btnId) => {
        let query = btnId.includes('#') ? btnId : '#' + sketcherName + btnId
        let btns = jquery.find(query)
        let btn = btns[0]
        btn.click()
        return 1
    }

    const chemdoodleEditorBtn = ($btn) => {
        var btnId = $btn.attr('href')
        chemdoodleClick(btnId)
        return false
    }

    const chemdoodleGetMol = (canvas) => {
        var mol = canvas.getMolecule()
        var molFile = ChemDoodle.writeMOL(mol)
        return molFile
    }

    const chemdoodleSetMol = (canvas, mol) => {
        molcd = ChemDoodle.readMOL(mol)
        canvas.loadMolecule(molcd)
        return false
    }

    const getEditorDimensions = () => {
        return [clientWidth, clientHeight]
    }

    const onWindowResize = () => {
        chemdoodleResize(sketcher, getEditorDimensions())
        chemdoodleClick('_button_center')
    }
</script>

<svelte:window on:resize={onWindowResize} />

<svelte:head>
    <script src="/chemdoodle/ChemDoodleWeb-unpacked.js"></script>
    <!-- <script src="/chemdoodle/ChemDoodleWeb-uis-unpacked.js" on:load={initializeChemdoodle}></script> -->
    <link rel="stylesheet" href="/chemdoodle/uis/jquery-ui-1.11.4.css" />
</svelte:head>

<div class="chemdoodle-editor p-5">
    <div class="shadow-outline rounded bg-white" style="width:500px; height:500px">
        <div class="chemdoodle-container h-full w-full" ref="chemdoodleCanvas" bind:clientWidth bind:clientHeight>
            <canvas id={sketcherName} />
            <div class="chemdoodle-hack1 block absolute bg-white " />
            <div class="chemdoodle-hack2 " />
        </div>

        <div class="chemdoodle-tools relative">
            <ul class="flex">
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_erase')} icon="eraser" tip="delete atom" />
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_undo')} icon="undo" tip="undo action" />
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_redo')} icon="redo" tip="redo action" />
                </li>
                <li>
                    <IconBtn
                        on:click={() => chemdoodleClick('_button_center')}
                        icon="compress-arrows-alt"
                        tip="re-center" />
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_scale_plus')} icon="search-plus" tip="zoom in" />
                </li>
                <li>
                    <IconBtn
                        on:click={() => chemdoodleClick('_button_scale_minus')}
                        icon="search-minus"
                        tip="zoom out" />
                </li>
            </ul>

            <br />

            <ul class="">
                {#each atoms as atom}
                    <IconBtn on:click={() => chemdoodleClick('_button_label_' + atom.toLowerCase())}>{atom}</IconBtn>
                {/each}
            </ul>

            <br />

            <ul class="">
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_bond_single')}>
                        <img
                            alt="single bond"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAMklEQVR42mNgGMSAm5qGOQDxfSDmopZhr4DYftSwUcNGDRs1jP6GcUPLM6oYBgNUKRwBiE8XjxDJvZUAAAAASUVORK5CYII=" />
                    </IconBtn>
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_bond_double')}>
                        <img
                            alt="double bond"
                            width="20"
                            height="20"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAQUlEQVR42mNgoA7gZqAicADi+0DMRS3DXgGxPTVcim4YzKVUMwyZP2rYqGGjhg0Dw7ihpQRVDIMBLmoahsulZAEA2GgvCVlTJIIAAAAASUVORK5CYII=" />
                    </IconBtn>
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_bond_triple')}>
                        <img
                            alt="Other Bond"
                            width="20"
                            height="20"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAUUlEQVR42mNgoAxwM1AROADxfSDmopZhr4DYnhouRTeMIpdiM4xslxIyDOZSqhmGzB81bNSwATeMG5rCqWIYDHBR0zBiXUoW4KKmYbhcShIAAA2MPiFy45L3AAAAAElFTkSuQmCC" />
                    </IconBtn>
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_ring_cyclobutane')}>
                        <img
                            alt="cyclobutane"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAiklEQVR42q3UwQ2AIAwF0E7hOv4bTCN7uwIXG0MTYgRa2iYNF3lB+ECkr4MCK3NX7hSBgfvmvtp4RmCCnB4Ug8lbKBaTTCiUH6tQGH9nimJzw3/R7IyEoEluQG0581RpzlvJuUL0KyRnaKd7b0VVB6lFTalYoVsRG6HwHOAXRcQTJmiJwPqcql7sB1sQMyMuYZLDAAAAAElFTkSuQmCC" />
                    </IconBtn>
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_ring_cyclohexane')}>
                        <img
                            alt="cyclohexane"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAoElEQVR42mNgIB4IMVAJSADxZCD+A8StQMxFiUH9QPwWiPuA2AiIlwLxYyCOAWJGcg2SQJO3AOITQHwSiC0pMQgZMEJd+Rjqall0BZxA/J8Ig9ABKDwboXo50SX/UxBx/4kWHDVw1EDcemEJu5+aCVsCmlPeEmEwwaxHisGW0ILhBLSgIKn4QjYYVHwtA+JHpBZf2AyeBC1gWygpYMmqAgB+TzRkG9cEtwAAAABJRU5ErkJggg==" />
                    </IconBtn>
                </li>
                <li>
                    <IconBtn on:click={() => chemdoodleClick('_button_ring_benzene')}>
                        <img
                            alt="benzene"
                            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAA5klEQVR42s2VPwrCMBSH06k0o1MX6R2c7A08hDfRRdDVgp5BR71GvUKdtHNvIMRf4FcoIdgkzeCDD0r+fLyQl1ch3GMmIkUOzuADDkBOEVWgA0ewABfQgjVIQkU5x/vMlqAGD1CGiFKwBS+QcSxhli2znpsyvVAZIh0r8AQ3UFiS0FnvuDczJ9XgW2++g4bSsVC/BiWPt+FxXSRqbEHmmZVyTvtvhTKmUMvesS+lYO0Fl01f2JWlsJvQws75UjpDnPL4Xk/PRdzLSjaGmo3Cq30Nxbp9XXlhXu3LJj6xwe6nNNigX8AXVupH9hGtsNcAAAAASUVORK5CYII=" />
                    </IconBtn>
                </li>
            </ul>
        </div>
    </div>
</div>

<style lang="postcss">
    /* chemdoodle  */

    .chemdoodle-container div:first-of-type {
        display: none;
    }

    .chemdoodle-container {
        position: relative;
    }

    .chemdoodle-hack1 {
        position: absolute;
        top: 10px;
        right: 10px;
        width: 20px;
        height: 20px;
        display: block;
        background: white;
    }

    .chemdoodle-tools li {
        @apply pl-2 pt-2;
   }
</style>
