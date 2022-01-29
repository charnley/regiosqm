<script>
    import Icon from './Icon.svelte'
</script>

<svelte:head />

<!-- This example requires Tailwind CSS v2.0+ -->
<div class="relative">
    <div class="max-w-7xl mx-auto px-4 sm:px-6">
        <div class="flex justify-between items-center border-b-2 border-gray-100 py-6 md:justify-start md:space-x-10">
            <div class="flex justify-start lg:w-0 lg:flex-1">
                <a class="text-2xl" href="#beginning">
                    <span class="space"><Icon type="fas" icon="atom" /></span>
                    <span class="" style="padding-left:5px;">RegioSQM</span>
                </a>
            </div>
            <nav class="hidden md:flex space-x-10">
                <!-- <a href="#" class="text-base font-medium text-gray-500 hover:text-gray-900"> -->
                <!--   Pricing -->
                <!-- </a> -->
                <!-- <a href="#" class="text-base font-medium text-gray-500 hover:text-gray-900"> -->
                <!--   Docs -->
                <!-- </a> -->
            </nav>
            <div class="hidden md:flex items-center justify-end md:flex-1 lg:w-0">
                <a href="#about" class="whitespace-nowrap text-base font-medium text-gray-500 hover:text-gray-900">
                    About
                </a>
                <a
                    href="https://github.com/jensengroup/regiosqm"
                    class="ml-8 whitespace-nowrap items-center justify-center px-4 py-2 border border-transparent rounded-md shadow-sm text-base font-medium text-white bg-indigo-600 hover:bg-indigo-700">
                    <Icon type="fab" icon="github" />
                    <span>GitHub</span>
                </a>
            </div>
        </div>
    </div>
</div>

<style>
</style>
