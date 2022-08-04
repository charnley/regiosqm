<script>
    import {Meta, Template, Story} from '@storybook/addon-svelte-csf'
    import IconBtn from './IconBtn.svelte'
</script>

<Meta
    title="Example/Button"
    component={IconBtn}
    argTypes={{
        label: {control: 'text'},
        tip: {control: 'text'},
        icon: {control: 'text'},
        primary: {control: 'boolean'},
        onClick: {action: 'onClick'},
    }} />

<Template let:args>
    <IconBtn tip={args.tip} icon={args.icon}>{args.label}</IconBtn>
</Template>

<Story
    name="Carbon Atom"
    args={{
        label: 'C',
    }} />

<Story
    name="Delete example"
    args={{
        label: "",
        icon: 'search-minus',
        tip: 'zoom out',
    }} />
