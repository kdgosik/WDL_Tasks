'use strict'

const IPFS = require('ipfs')

const ipfs = new IPFS({
  EXPERIMENTAL: {
    pubsub: true
  }
})

// example obj
const obj = {
  a: 1,
  b: [1, 2, 3],
  c: {
    ca: [5, 6, 7],
    cb: 'foo'
  }
}

ipfs.dag.put(obj, { format: 'dag-cbor', hashAlg: 'sha2-256' }, (err, cid) => {
  console.log(cid.toBaseEncodedString())
  // zdpuAmtur968yprkhG9N5Zxn6MFVoqAWBbhUAkNLJs2UtkTq5
})

function errOrLog(err, result) {
  if (err) {
    console.error('error: ' + err)
  } else {
    console.log(result)
  }
}

ipfs.dag.tree('zdpuAmtur968yprkhG9N5Zxn6MFVoqAWBbhUAkNLJs2UtkTq5', errOrLog)
