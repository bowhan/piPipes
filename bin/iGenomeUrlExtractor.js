// this code parse the igenome page and extract the hrefs ending with .tar.gz and the its name
// some of the names are not ready to be used as a genome name by piPipes, need manual inspection
// Bo Han 2016-09-01
// run it as node iGenomeUrlExrtractor.js after installing jsdom (npm install jsdom)

var jsdom = require("jsdom");

jsdom.env({
              url: 'http://support.illumina.com/sequencing/sequencing_software/igenome.html',
              done: function (err, window) {
                  global.window   = window;
                  global.document = window.document;
                  extractHrefs();
              }
          });

function extractHrefs() {
    let tdas = document.querySelectorAll('td > a');
    for (let i = 0; i < tdas.length; ++i) {
        let x = tdas[i];
        if (x.href.endsWith('.tar.gz'))
            console.log(`${x.text}=${x.href}`);
    }
}