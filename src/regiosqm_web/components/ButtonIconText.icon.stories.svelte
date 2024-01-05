<script>
    import {Meta, Template, Story} from '@storybook/addon-svelte-csf'
    import ButtonIconText from './ButtonIconText.svelte'
    import Tile from 'carbon-components-svelte/src/Tile/Tile.svelte'
    import Button from 'carbon-components-svelte/src/Button/Button.svelte'
</script>

<Meta
    title="CarbonButton IconSlot"
    component={ButtonIconText}
    argTypes={{
        iconDescription: {control: 'text'},
        kind: {control: 'text'},
    }} />

<Template let:args>
    <Tile>
        <ButtonIconText iconDescription={args.iconDescription} kind={args.kind}>
            <img
                alt="benzene"
                src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABQAAAAUCAYAAACNiR0NAAAA5klEQVR42s2VPwrCMBSH06k0o1MX6R2c7A08hDfRRdDVgp5BR71GvUKdtHNvIMRf4FcoIdgkzeCDD0r+fLyQl1ch3GMmIkUOzuADDkBOEVWgA0ewABfQgjVIQkU5x/vMlqAGD1CGiFKwBS+QcSxhli2znpsyvVAZIh0r8AQ3UFiS0FnvuDczJ9XgW2++g4bSsVC/BiWPt+FxXSRqbEHmmZVyTvtvhTKmUMvesS+lYO0Fl01f2JWlsJvQws75UjpDnPL4Xk/PRdzLSjaGmo3Cq30Nxbp9XXlhXu3LJj6xwe6nNNigX8AXVupH9hGtsNcAAAAASUVORK5CYII=" />
        </ButtonIconText>
    </Tile>
</Template>

<Story
    name="Benzene"
    args={{
        iconDescription: 'Benzene',
    }} />
