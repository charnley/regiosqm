<script context="module">
    import Icon from './Icon.svelte'
</script>

<script>
    export let icon = undefined
    export let tip = undefined
</script>

<!-- svelte-ignore a11y-missing-attribute -->
<a on:click title={tip}><Icon {icon} /><slot class="content" /></a>

<style>
    a {
        display: block;
        width: 20px;
        height: 20px;
        background-color: grey;
    }
</style>
