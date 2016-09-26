/*
	# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
	# Copyright (C) 2014-2016  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
	#
	# This program is free software; you can redistribute it and/or modify
	# it under the terms of the GNU General Public License as published by
	# the Free Software Foundation; either version 3 of the License, or
	# (at your option) any later version.

	# This program is distributed in the hope that it will be useful,
	# but WITHOUT ANY WARRANTY; without even the implied warranty of
	# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	# GNU General Public License for more details.

	# You should have received a copy of the GNU General Public License along
	# with this program; if not, write to the Free Software Foundation, Inc.,
	# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    This javascript/node code extracts the URL of iGenome files and assign them to a variable in bash format
 
* usage:
*   npm install jsdom
*   node iGenomeUrlExtractor.js
* */

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
    let tdas = document.querySelectorAll('td a');
    for (let i = 0; i < tdas.length; ++i) {
        let x = tdas[i];
        if (x.href.endsWith('.tar.gz'))
            console.log(`${x.text}=${x.href}`);
    }
}