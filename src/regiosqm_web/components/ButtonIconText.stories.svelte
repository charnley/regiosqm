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
        text: {control: 'text'},
        kind: {control: 'text'},
    }} />

<Template let:args>
    <Tile>
        <ButtonIconText iconDescription={args.iconDescription} kind={args.kind}>{args.text}</ButtonIconText>
    </Tile>
</Template>

<Story
    name="Carbon Atom"
    args={{
        text: 'C',
        iconDescription: 'Add Carbon',
    }} />

<Story
    name="Charge"
    args={{
        text: '+1',
        iconDescription: 'Add Charge',
    }} />
