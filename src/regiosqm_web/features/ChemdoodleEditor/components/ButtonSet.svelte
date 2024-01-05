<script>
    // your script goes here
</script>

<div class="chemdoodle-button-set">
    <slot />
</div>

<style>
    .chemdoodle-button-set {
        background: red;
    }
    :global(.chemdoodle-button-set button) {
        margin: 0 5px 5px 0;
    }
</style>
