<script>
    export let icon
    export let type = 'fas'

    let tag = 'fa-' + icon
</script>

{#if tag}
    <i class="{type} {tag}" />
{/if}

<style>
</style>
