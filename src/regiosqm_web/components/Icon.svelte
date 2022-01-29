<script>
    export let icon

    let tag = 'fa-' + icon
</script>

{#if tag}
    <i class="fas {tag}" />
{/if}

<style>
</style>
