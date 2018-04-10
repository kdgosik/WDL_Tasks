'use strict'

var ipfsAPI = require('ipfs-api')

// connect to ipfs daemon API server
var ipfs = ipfsAPI('localhost', '5001', {protocol: 'http'}) // leaving out the arguments will default to these values


ipfs.id(function (err, identity) {
  if (err) {
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
  }

  ipfs.id(function(err, ident) {
    if (err) {
      throw err
    }
    console.log(ident)
  })
})
