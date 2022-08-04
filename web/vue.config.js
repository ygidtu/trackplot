

module.exports = {
    lintOnSave: true, //禁用eslint
    publicPath: 'static',  // static
    outputDir: "../ui",
    runtimeCompiler: true,
    productionSourceMap: false,
    parallel: 4,
    configureWebpack:{
        performance: {
            hints: false
        },
        optimization: {
            splitChunks: {
                maxSize: 500000,
            }
        }
    },
    // options...
    devServer: {
        disableHostCheck: true,
    },
};