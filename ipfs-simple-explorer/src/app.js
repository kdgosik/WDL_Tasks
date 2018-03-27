'use strict'

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

document.addEventListener('DOMContentLoaded', () => {
  document.getElementById('store').onclick = store
})
