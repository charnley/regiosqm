import '../public/build/bundle.css'
import '../public/fonts/font.css'
import '../public/fontawesome/css/all.css'

export const parameters = {
    actions: {argTypesRegex: '^on[A-Z].*'},
    controls: {
        matchers: {
            color: /(background|color)$/i,
            date: /Date$/,
        },
    },
}
