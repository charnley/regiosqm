<script context="module">
    let sketcher

    export const chemdoodleGetMol = () => {
        let mol = sketcher.getMolecule()
        let molFile = ChemDoodle.writeMOL(mol)
        return molFile
    }

    export const chemdoodleSetMol = (mol) => {
        let molcd = ChemDoodle.readMOL(mol)
        sketcher.loadMolecule(molcd)
    }
</script>

<script>
    import {jquery, chemdoodle} from '../../stores/libs.js'

    import TrashCan from 'carbon-icons-svelte/lib/TrashCan.svelte'
    import Undo from 'carbon-icons-svelte/lib/Undo.svelte'
    import Redo from 'carbon-icons-svelte/lib/Redo.svelte'
    import CenterToFit from 'carbon-icons-svelte/lib/CenterToFit.svelte'
    import ZoomIn from 'carbon-icons-svelte/lib/ZoomIn.svelte'
    import ZoomOut from 'carbon-icons-svelte/lib/ZoomOut.svelte'

    import {Grid, Row, Column} from 'carbon-components-svelte'
    import {Button} from 'carbon-components-svelte'

    import ButtonIconText from '../../components/ButtonIconText.svelte'
    import ButtonSet from './components/ButtonSet.svelte'

    export let sketcherName = 'sketcherSingle'
    export let atoms = ['H', 'Li', 'Be', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Br', 'I']

    let clientWidth
    let clientHeight

    const initializeChemdoodle = () => {
        sketcher = new $chemdoodle.SketcherCanvas(sketcherName, 100, 100, {useServices: false, oneMolecule: true})
        chemdoodleResize(sketcher, [500, 500])
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
        let btns = $jquery.find(query)
        let btn = btns[0]
        btn.click()
    }

    const getEditorDimensions = () => {
        return [clientWidth, clientHeight]
    }

    const onWindowResize = () => {
        chemdoodleResize(sketcher, getEditorDimensions())
        chemdoodleClick('_button_center')
    }

    // Wait for chemdoodle
    $: $chemdoodle && initializeChemdoodle()

    const btnOptions = {size: 'sm'}
</script>

<svelte:window on:resize={onWindowResize} />

<Grid fullWidth narrow>
    <Row>
        <Column>
            <Button
                {...btnOptions}
                on:click={() => chemdoodleClick('_button_erase')}
                icon={TrashCan}
                iconDescription="Delete Atom" />
            <Button
                {...btnOptions}
                on:click={() => chemdoodleClick('_button_undo')}
                icon={Undo}
                iconDescription="Undo last action" />
            <Button
                {...btnOptions}
                on:click={() => chemdoodleClick('_button_redo')}
                icon={Redo}
                iconDescription="Redo last undo" />
            <Button
                {...btnOptions}
                on:click={() => chemdoodleClick('_button_center')}
                icon={CenterToFit}
                iconDescription="Re-center editor" />
            <Button
                size="sm"
                on:click={() => chemdoodleClick('_button_scale_plus')}
                icon={ZoomIn}
                iconDescription="Zoom in" />
            <Button
                size="sm"
                on:click={() => chemdoodleClick('_button_scale_minus')}
                icon={ZoomOut}
                iconDescription="Zoom out" />

            <div>
                <ButtonSet>
                    {#each atoms as atom}
                        <ButtonIconText
                            {...btnOptions}
                            size="sm"
                            on:click={() => chemdoodleClick('_button_label_' + atom.toLowerCase())}
                            >{atom}</ButtonIconText>
                    {/each}
                </ButtonSet>
            </div>

            <div>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_bond_single')}>
                    <img
                        alt="single bond"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAMklEQVR42mNgGMSAm5qGOQDxfSDmopZhr4DYftSwUcNGDRs1jP6GcUPLM6oYBgNUKRwBiE8XjxDJvZUAAAAASUVORK5CYII=" />
                </ButtonIconText>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_bond_double')}>
                    <img
                        alt="double bond"
                        width="16"
                        height="16"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAQUlEQVR42mNgoA7gZqAicADi+0DMRS3DXgGxPTVcim4YzKVUMwyZP2rYqGGjhg0Dw7ihpQRVDIMBLmoahsulZAEA2GgvCVlTJIIAAAAASUVORK5CYII=" />
                </ButtonIconText>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_bond_triple')}>
                    <img
                        alt="Other Bond"
                        width="16"
                        height="16"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAUUlEQVR42mNgoAxwM1AROADxfSDmopZhr4DYnhouRTeMIpdiM4xslxIyDOZSqhmGzB81bNSwATeMG5rCqWIYDHBR0zBiXUoW4KKmYbhcShIAAA2MPiFy45L3AAAAAElFTkSuQmCC" />
                </ButtonIconText>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_ring_cyclobutane')}>
                    <img
                        width="16"
                        height="16"
                        alt="cyclobutane"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAiklEQVR42q3UwQ2AIAwF0E7hOv4bTCN7uwIXG0MTYgRa2iYNF3lB+ECkr4MCK3NX7hSBgfvmvtp4RmCCnB4Ug8lbKBaTTCiUH6tQGH9nimJzw3/R7IyEoEluQG0581RpzlvJuUL0KyRnaKd7b0VVB6lFTalYoVsRG6HwHOAXRcQTJmiJwPqcql7sB1sQMyMuYZLDAAAAAElFTkSuQmCC" />
                </ButtonIconText>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_ring_cyclohexane')}>
                    <img
                        width="16"
                        height="16"
                        alt="cyclohexane"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAAoElEQVR42mNgIB4IMVAJSADxZCD+A8StQMxFiUH9QPwWiPuA2AiIlwLxYyCOAWJGcg2SQJO3AOITQHwSiC0pMQgZMEJd+Rjqall0BZxA/J8Ig9ABKDwboXo50SX/UxBx/4kWHDVw1EDcemEJu5+aCVsCmlPeEmEwwaxHisGW0ILhBLSgIKn4QjYYVHwtA+JHpBZf2AyeBC1gWygpYMmqAgB+TzRkG9cEtwAAAABJRU5ErkJggg==" />
                </ButtonIconText>
                <ButtonIconText {...btnOptions} on:click={() => chemdoodleClick('_button_ring_benzene')}>
                    <img
                        width="16"
                        height="16"
                        alt="benzene"
                        src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAA5klEQVR42s2VPwrCMBSH06k0o1MX6R2c7A08hDfRRdDVgp5BR71GvUKdtHNvIMRf4FcoIdgkzeCDD0r+fLyQl1ch3GMmIkUOzuADDkBOEVWgA0ewABfQgjVIQkU5x/vMlqAGD1CGiFKwBS+QcSxhli2znpsyvVAZIh0r8AQ3UFiS0FnvuDczJ9XgW2++g4bSsVC/BiWPt+FxXSRqbEHmmZVyTvtvhTKmUMvesS+lYO0Fl01f2JWlsJvQws75UjpDnPL4Xk/PRdzLSjaGmo3Cq30Nxbp9XXlhXu3LJj6xwe6nNNigX8AXVupH9hGtsNcAAAAASUVORK5CYII=" />
                </ButtonIconText>
            </div>
        </Column>
        <Column>
            <div class="chemdoodle-container h-full w-full" ref="chemdoodleCanvas" bind:clientWidth bind:clientHeight>
                <canvas id={sketcherName} />
                <div class="chemdoodle-hack1 block absolute bg-white " />
                <div class="chemdoodle-hack2 " />
            </div>
        </Column>
    </Row>
</Grid>

<!-- <div class="chemdoodle-editor p-5"> -->
<!--     <div class="shadow-outline rounded bg-white relative overflow-hidden" style=""> -->
<!--         <div class=" h-96 " /> -->

<!--         <div class="chemdoodle-tools relative"> -->
<!--             <ul class=""> -->
<!--                 <li> -->
<!--                 </li> -->
<!--                 <li> -->
<!--                 </li> -->
<!--                 <li> -->
<!--                 </li> -->
<!--                 <li> -->
<!--                 </li> -->
<!--                 <li> -->
<!--                 </li> -->
<!--                 <li> -->
<!--                 </li> -->
<!--             </ul> -->

<!--             <ul class=""> -->
<!--                 {#each atoms as atom} -->
<!--                     <li> -->
<!--                         <ButtonIconText on:click={() => chemdoodleClick('_button_label_' + atom.toLowerCase())} -->
<!--                             >{atom}</ButtonIconText> -->
<!--                     </li> -->
<!--                 {/each} -->
<!--             </ul> -->

<!--             <ul class=""> -->
<!--                 <li> -->
<!--                 </li> -->
<!--             </ul> -->
<!--         </div> -->
<!--     </div> -->

<!-- </div> -->
<style lang="postcss">
    /* chemdoodle  */

    :global(.chemdoodle-container div:first-of-type) {
        display: none;
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

    .chemdoodle-container {
        @apply relative bg-gray-400;
    }

    .chemdoodle-tools {
        @apply relative bg-slate-100;
    }
    .chemdoodle-tools ul {
        @apply flex p-2;
    }

    .chemdoodle-tools li {
        @apply pl-2;
    }
</style>
