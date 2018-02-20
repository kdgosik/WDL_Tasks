const graphql = require('graphql');
const express = require('express');
const expressGraphQL = require('express-graphql');
const schema = require('./schema/schema_server.js');

const app = express();

app.use('/graphql', expressGraphQL({
  schema:schema,
  graphiql:true
}));

app.listen(4000, () => {
  console.log("Server is running on port 4000..")
});
