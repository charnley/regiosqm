<script context="module">
    import Button from 'carbon-components-svelte/src/Button/Button.svelte'
</script>

<script>
    export let iconDescription = undefined

    $: hackClass = iconDescription
        ? 'bx--tooltip__trigger bx--tooltip--a11y bx--tooltip--align-center bx--btn--icon-only--bottom'
        : ''
    const standardClass = 'bx--btn--icon-only bx--btn--icontext'
</script>

<Button class="{hackClass} {standardClass}" {...$$restProps}>
    <span><slot /> </span>
    {#if iconDescription}
        <span class="bx--assistive-text">{iconDescription}</span>
    {/if}
</Button>

<style>
    :global(.bx--btn--icontext span) {
        display: block;
        height: 16px;
        width: 16px;
        text-align: center;
    }

    :global(.bx--btn--icontext img) {
        display: inline;
        filter: invert(100%);
        height: 16px;
        line-height: 16px;
        width: 16px;
    }
</style>
