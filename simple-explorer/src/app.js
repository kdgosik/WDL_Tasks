'use strict'
/*
// this part cretes its own ifps intance

const IPFS = require('ipfs')

const ipfs = new IPFS({
  EXPERIMENTAL: {
    pubsub: true
  }
})


ipfs.once('ready', () => ipfs.id((err, info) => {
  if (err) { throw err }
  console.log('IPFS node ready with address ' + info.id)
}))
*/

var ipfsAPI = require('ipfs-api')

// connect to ipfs daemon API server
var ipfs = ipfsAPI('localhost', '5001', {protocol: 'http'}) // leaving out the arguments will default to these values


function store () {
  var toStore = document.getElementById('source').value

  ipfs.files.add(Buffer.from(toStore), (err, res) => {
    if (err || !res) {
      return console.error('ipfs add error', err, res)
    }

    res.forEach((file) => {
      if (file && file.hash) {
        console.log('successfully stored', file.hash)
        display(file.hash)
      }
    })
  })
}


function display (hash) {
  // buffer: true results in the returned result being a buffer rather than a stream
  ipfs.files.cat(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash').innerText = hash
    document.getElementById('content').innerText = data
  })
}

// TO BE REMOVED

function store1 () {
  var toStore = document.getElementById('source1').value

  ipfs.files.add(Buffer.from(toStore), (err, res) => {
    if (err || !res) {
      return console.error('ipfs add error', err, res)
    }

    res.forEach((file) => {
      if (file && file.hash) {
        console.log('successfully stored', file.hash)
        display1(file.hash)
      }
    })
  })
}


function display1 (hash) {
  // buffer: true results in the returned result being a buffer rather than a stream
  ipfs.files.cat(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash1').innerText = hash
    document.getElementById('content1').innerText = data
  })
}


// TO BE REMOVED
// needs lots of workflow
// https://github.com/ipfs/interface-ipfs-core/blob/master/SPEC/OBJECT.md#objectpatchaddlink
function store2 () {
  var toStore = document.getElementById('source').value
  var hash1 = document.getElementById('hash1').innerText
  console.log(hash1)


  ipfs.files.add(Buffer.from(toStore), (err, node1) => {
    if (err || !node1) {
      return console.error('ipfs add error', err, node1)
    }
  console.log(node1[0].hash)
    ipfs.object.patch.addLink(node1[0].hash, {
      name:'link-to',
      size: 10,
      multihash: hash1
    }, (err, res) => {
      if (err || !res) {
        return console.error('ipfs object error', err, res)
      }
      console.log(res.multihash)
      display2(res.multihash)

      /*
      res.forEach((file) => {
        if (file && file.hash) {
          console.log('successfully stored', file.hash)
          display2(file.hash)
        }
      })
      */

    })
  })
}


function display2 (hash) {
  ipfs.object.get(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash2').innerText = data.toJSON().multihash
    document.getElementById('content2').innerText = JSON.stringify(data.toJSON()) + '\n\n' + data.toJSON().data
  })
}


document.addEventListener('DOMContentLoaded', () => {
  document.getElementById('store').onclick = store
  document.getElementById('store1').onclick = store1
  document.getElementById('store2').onclick = store2
})


function httpGetAsync(theUrl, callback) {
    var xmlHttp = new XMLHttpRequest();
    xmlHttp.onreadystatechange = function() {
        if (xmlHttp.readyState == 4 && xmlHttp.status == 200)
            callback(JSON.parse(this.responseText));
    }


    xmlHttp.open("GET", theUrl, true); // true for asynchronous
    xmlHttp.setRequestHeader("Content-type", "application/json");
    xmlHttp.send(null);
}

// make function to fill in wdl list select

/*
Golang Code
content, err := ioutil.ReadFile("../wdl_hash.json")
if err != nil {
  fmt.Println(err)
}
'https://api.firecloud.org/ga4gh/v1/tools'
*/

function loadSelect() {

  var id = [];
  var toolname = [];
  var version = [];

  httpGetAsync('https://api.firecloud.org/ga4gh/v1/tools', function(result) {
    for (let i=0; i<result.length; i++) {
      id.push(result[i].id)
      toolname.push(result[i].toolname)
      version.push(result[i]["meta-version"])
    }

    var x1 = document.getElementById("wdl1");
    for (let i=0; i<result.length; i++) {
      var c = document.createElement("option");
      c.text = toolname[i];
      c.value =id[i];
      x1.options.add(c);
    }

    var x2 = document.getElementById("wdl2");
    for (let i=0; i<result.length; i++) {
      var c = document.createElement("option");
      c.text = toolname[i];
      c.value =id[i];
      x2.options.add(c);
    }
  })
}

// use above function to load selections
loadSelect()




// Listen for form submit
document.getElementById('wdl_data1').addEventListener('submit', submitForm);
// Listen for form submit
document.getElementById('wdl_data2').addEventListener('submit', submitForm);

// Submit form
function submitForm(e){
  e.preventDefault();

  // Get values
  httpGetAsync("https://api.firecloud.org/ga4gh/v1/tools/"+e.value+"/versions/"+val.MetaVersion+"/plain-WDL/descriptor", function(result) {

  // Show alert
  document.querySelector('.alert').style.display = 'block';

  // Hide alert after 3 seconds
  setTimeout(function(){
    document.querySelector('.alert').style.display = 'none';
  },3000);

  // Clear form
  document.getElementById('contactForm').reset();
}
