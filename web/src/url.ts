export function createURL(...urls: string[]) {
    let res = [];

    // res.push("http://127.0.0.1:5000")

    for (let i=0; i < urls.length; i++) {
        let temp = urls[i].toString();
        temp = temp.replace(/(^\/|\/$)/g, "");

        if (temp !== "") {
            res.push(temp)
        }
    }
    return res.join("/")
}

const urls = {
    static: createURL("/static"),
    del: createURL('/api/del'),
    file: createURL('/api/file'),
    params: createURL('/api/params'),
    plot: createURL('/api/plot'),
    log: createURL('/api/log')
};

export default urls;
