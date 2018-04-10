'use strict'

// this part creates its own ifps intance

// make conditional statement to ping localhost:5001 for the API
// if error then start one in the browser


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

/*
var ipfsAPI = require('ipfs-api')

// connect to ipfs daemon API server
var ipfs = ipfsAPI('localhost', '5001', {protocol: 'http'}) // leaving out the arguments will default to these values
*/


function store1 () {
  var obj = {};
  var id = document.getElementById("wdl1").value
  var toStore = document.getElementById('source1').value

  obj["id"] = id;
  obj["data"] = toStore;
  console.log(obj.id)
  console.log(obj.data)
  console.log(JSON.stringify(obj))

// add object to ipfs
  ipfs.files.add(Buffer.from(JSON.stringify(obj)), (err, res) => {
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
  ipfs.object.get(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash1').innerText = hash
    // document.getElementById('source1').innerText = data.toJSON().data
  })
}


function store2 () {
  var obj = {};
  var id = document.getElementById("wdl2").value
  var toStore = document.getElementById('source2').value

  obj["id"] = id;
  obj["data"] = toStore;
  console.log(obj.id)
  console.log(obj.data)
  console.log(JSON.stringify(obj))

// add object to ipfs
  ipfs.files.add(Buffer.from(JSON.stringify(obj)), (err, res) => {
    if (err || !res) {
      return console.error('ipfs add error', err, res)
    }

    res.forEach((file) => {
      if (file && file.hash) {
        console.log('successfully stored', file.hash)
        display2(file.hash)
      }
    })
  })
}


function display2 (hash) {
  ipfs.object.get(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash2').innerText = hash
    // document.getElementById('source2').innerText = data.toJSON().data
  })
}



function store3 () {
  var hash1 = document.getElementById('hash1').value
  console.log(hash1)
  var hash2 = document.getElementById('hash2').value
  console.log(hash2)

  ipfs.object.patch.addLink(hash1, {
      name:'link-to',
      multihash: hash2,
      size: 10
    }, (err, res) => {
      if (err || !res) {
        return console.error('ipfs object error', err, res)
      }
      console.log(res.toJSON())
      console.log(res.toJSON().multihash)
      console.log(res.multihash)
      display3(res.multihash)
    })
}


function display3 (hash) {
  ipfs.object.get(hash, (err, data) => {
    if (err) { return console.error('ipfs cat error', err) }

    document.getElementById('hash3').value = data.toJSON().multihash
    document.getElementById('hash3').innerText = data.toJSON().multihash
  //  document.getElementById('content3').innerText = JSON.stringify(data.toJSON()) + '\n\n' + data.toJSON().data
  })
}


document.addEventListener('DOMContentLoaded', () => {
  document.getElementById('store1').onclick = store1
  document.getElementById('store2').onclick = store2
  document.getElementById('store3').onclick = store3
})


function httpGetAsync(theUrl, headerType, callback) {
    var xmlHttp = new XMLHttpRequest();
    xmlHttp.onreadystatechange = function() {
        if (xmlHttp.readyState == 4 && xmlHttp.status == 200)
          if(headerType == "application/json"){
            callback(JSON.parse(this.responseText));
          }
          else{
            callback(this.responseText);
          }
    }


    xmlHttp.open("GET", theUrl, true); // true for asynchronous
    xmlHttp.setRequestHeader("Content-type", headerType);
    xmlHttp.send(null);
}

// make function to fill in wdl list select
function loadSelect() {

  var id = [];
  var toolname = [];
  var version = [];

  httpGetAsync('https://api.firecloud.org/ga4gh/v1/tools', "application/json", function(result) {
    //var obj = JSON.parse(result)
    for (let i=0; i<result.length; i++) {
      id.push(result[i].id)
      toolname.push(result[i].toolname)
      version.push(result[i]["meta-version"])
    }

    var x1 = document.getElementById("wdl1");
    for (let i=0; i<result.length; i++) {
      var c = document.createElement("option");
      c.text = toolname[i];
      c.value =id[i]+"-"+version[i];
      x1.options.add(c);
    }

    var x2 = document.getElementById("wdl2");
    for (let i=0; i<result.length; i++) {
      var c = document.createElement("option");
      c.text = toolname[i];
      c.value =id[i]+"-"+version[i];
      x2.options.add(c);
    }
  })
}

// use above function to load selections
loadSelect()


// Listen for form submit
document.getElementById('wdl_data1').addEventListener('submit', submitForm1);
// Listen for form submit
document.getElementById('wdl_data2').addEventListener('submit', submitForm2);
// Listen for form submit
document.getElementById('display_link').addEventListener('submit', submitForm3);



// Submit form wdl_data1
function submitForm1(e){
  e.preventDefault();

  var formVal = document.getElementById("wdl1").value
  console.log(formVal)

  var toolId = formVal.split("-")[0];
  console.log(toolId)
  var toolVersion = formVal.split("-")[1];
  console.log(toolVersion)

  // Get values
  httpGetAsync("https://api.firecloud.org/ga4gh/v1/tools/"+toolId+"/versions/"+toolVersion+"/plain-WDL/descriptor", "text/plain", function(result) {

    document.getElementById('source1').innerText = result

  })
  // Show alert
  document.querySelector('.alert').style.display = 'block';

  // Hide alert after 3 seconds
  setTimeout(function(){
    document.querySelector('.alert').style.display = 'none';
  },3000);

  // Clear form
  //document.getElementById('contactForm').reset();
}


// Submit form wdl_data2
function submitForm2(e){
  e.preventDefault();

  var formVal = document.getElementById("wdl2").value
  console.log(formVal)

  var toolId = formVal.split("-")[0];
  console.log(toolId)
  var toolVersion = formVal.split("-")[1];
  console.log(toolVersion)

  // Get values
  httpGetAsync("https://api.firecloud.org/ga4gh/v1/tools/"+toolId+"/versions/"+toolVersion+"/plain-WDL/descriptor", "text/plain", function(result) {

    document.getElementById('source2').innerText = result

  })
  // Show alert
  document.querySelector('.alert').style.display = 'block';

  // Hide alert after 3 seconds
  setTimeout(function(){
    document.querySelector('.alert').style.display = 'none';
  },3000);

  // Clear form
  //document.getElementById('contactForm').reset();
}



// Submit form display_link
function submitForm3(e){
  e.preventDefault();

  var linkhash = document.getElementById("linkhash").value;
  console.log(linkhash)

  ipfs.object.get(linkhash, (err, result1) => {
    if (err || !result1) {
      return console.error('ipfs add error', err, result1)
    }

    document.getElementById('outwdl1').innerText = result1.toJSON().data
    console.log(result1.toJSON().multihash)
    console.log(result1.toJSON().links)

    result1.toJSON().links.forEach((result2) => {
      ipfs.object.get(result2.multihash, (err, result3) => {
        if (err || !result3) {
          return console.error('ipfs add error', err, result3)
        }

        document.getElementById('outwdl2').innerText = result3.toJSON().data
      })
    })



  })

  // Show alert
  document.querySelector('.alert').style.display = 'block';

  // Hide alert after 3 seconds
  setTimeout(function(){
    document.querySelector('.alert').style.display = 'none';
  },3000);

  // Clear form
  //document.getElementById('contactForm').reset();
}
