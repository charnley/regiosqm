<script context="module">
    import Icon from './Icon.svelte'
</script>

<script>
    export let icon = undefined
    export let tip = undefined
</script>

<!-- svelte-ignore a11y-missing-attribute -->
<a on:click title={tip}><Icon {icon} /><slot class="content" /></a>

<style lang="postcss">
    a {
        height: 30px;
        width: 30px;
        line-height: 30px;
        padding: 0;
        text-align: center;
        font-size: 15px;
        font-weight: 700;

        user-select: none;
        cursor: pointer;
        display: inline-block;
        border: 1px solid #ddd;
        border-color: #fff;
        border-radius: 3px;
        text-decoration: none;
        color: #222;
        background-color: white;
    }

    a:hover {
        @apply bg-blue-700;
    }

    :global(a:hover img) {
        filter: invert(100%);
    }
</style>
