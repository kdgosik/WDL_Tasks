# WDL_Tasks

## install

Install the necessary node modules

```
npm install --save
```

## Run Server

We will need two terminals to run this.  There are two parts that need ports to run.  The first is a simple json server that creates a json database on localhost:3000.  To do this by the json-server package. The package.json file has a script to do this by typing the following command.

```
npm run json:server
```

In a different terminal, we will start the graphiql server.  This will run on localhost:4000/graphql.  This will have the IDE provided called Graphiql.  The package.json file has a script to run this as well by using the following command.

```
npm run dev:server
```

Navigate to <a href="http://localhost:4000/graphql">http://localhost:4000/graphql</a> where you can execute queries.  Here is an example one to get started.

```
{
  task(id: 3){
    id
    name
    input
    output
    connections{
      id
      name
    }
  }
}
```
