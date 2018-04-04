var ipfsAPI = require('ipfs-api')

// connect to ipfs daemon API server
var ipfs = ipfsAPI('localhost', '5001', {protocol: 'http'}) // leaving out the arguments will default to these values

ipfs.files.add(Buffer.from('some data to store'), (err, res) => {
  if (err || !res) {
    return console.error('ipfs add error', err, res)
  }

  res.forEach((file) => {
    if (file && file.hash) {
      console.log('successfully stored', file.hash)
    }
  })
})
