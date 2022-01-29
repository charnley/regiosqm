<script>
    import Editor from './components/Editor.svelte'
    import Navigation from './components/Navigation.svelte'
    export let appTitle = 'Regioselectitiy App'
</script>

<svelte:head>
    <title>{appTitle}</title>
    <meta name="theme-color" content="#ffffff" />
    <meta name="google" content="notranslate" />
    <meta name="mobile-web-app-capable" content="yes" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0, maximum-scale=1.0, user-scalable=no" />
    <link href="https://fonts.googleapis.com/css?family=Roboto|Fira+Sans:300,400,500" rel="stylesheet" />
    <link href="https://fonts.googleapis.com/css?family=Nunito:400,700,800" rel="stylesheet" />
    <link rel="stylesheet" href="/fontawesome/css/all.css" />
</svelte:head>

<div class="min-h-screen">
    <Navigation />

    <section class="relative mt-20">
        <div class="container flex sm:flex-nowrap flex-wrap max-w-7xl mx-auto px-4 sm:px-6">
            <div class="lg:w-1/3 md:w-1/2 bg-white flex flex-col md:ml-auto w-full">
                <h1 class="text-4xl mb-1">Predict Regioselectitiy</h1>
                <p>
                    of electrophilic aromatic substitution reactions. For citation and method details, please see the
                    RegioSQM/ML papers.
                </p>
                <p><a href="https://fu">Read more...</a></p>
            </div>
            <div class="lg:w-2/3 md:w-1/2 bg-gray-300 rounded-lg sm:mr-10 relative">
                <Editor />
            </div>
        </div>
    </section>
</div>

<style global lang="postcss">
    @tailwind base;
    @tailwind components;
    @tailwind utilities;

    body {
        font-family: Nunito, sans-serif;
        font-weight: bold;
    }
</style>
