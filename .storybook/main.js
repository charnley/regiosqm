module.exports = {
    stories: ['../src/**/*.stories.mdx', '../src/**/*.stories.@(js|jsx|ts|tsx|svelte)'],
    staticDirs: ['../public'],
    addons: [
        '@storybook/addon-links',
        '@storybook/addon-essentials',
        '@storybook/addon-svelte-csf',
        '@storybook/addon-interactions',
    ],
    features: {
        postcss: true,
    },
    framework: '@storybook/svelte',
    core: {
        builder: '@storybook/builder-webpack4',
    },
}
