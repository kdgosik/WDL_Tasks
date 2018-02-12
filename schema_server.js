const axios = require('axios');
const {
  GraphQLObjectType,
  GraphQLString,
  GraphQLInt,
  GraphQLSchema,
  GraphQLList,
  GraphQLNonNull
} = require('graphql');

// Hardcoded Data
// const customers = [
//   {id:'1', name:'John Doe', email: 'jdoe@email.com', age: 35},
//   {id:'2', name:'Steve Smith', email: 'ssmith@email.com', age: 25},
//   {id:'3', name:'Sarah Williams', email: 'swilliams@email.com', age: 30},
// ];

// CustomerType
const TaskType = new GraphQLObjectType({
  name:'Task',
  fields:() => ({
    id: {type: GraphQLInt},
    name: {type: GraphQLString},
    input: {type: GraphQLList(GraphQLString)},
    command: {type: GraphQLString},
    output: {type: GraphQLString},
    connections: {
      type: GraphQLList(TaskType),
      description: "Finding related tasks",
      resolve(args){
        let out = [];
        let conn = axios.get('http://localhost:3000/connections')
        .then(res => res.data)
        .then(function(result){
          for(let i = 0; i < result.length; i++){
            if(result[i].from == args.id){
              axios.get('http://localhost:3000/tasks/'+result[i].to)
              .then(res => res.data)
              .then(function(o){
                out.push(o)
                console.log(out)
              });
            }
          }
          console.log(conn) // Promise{ <pending> }
          console.log(out)
      })

    }
  }
  })
});

// Root Query
const RootQuery = new GraphQLObjectType({
  name: 'RootQueryType',
  fields: {
    task: {
      type: TaskType,
      args: {
        id:{type: GraphQLInt}
      },
      resolve(parentValue, args){
        return axios.get('http://localhost:3000/tasks/'+args.id)
          .then(res => res.data);
      }
    },
    tasks: {
      type: new GraphQLList(TaskType),
      resolve(parentValue, args){
        //return customers;
        return axios.get('http://localhost:3000/tasks/')
          .then(res => res.data);
      }
    }
  }
});

// Mutations
const mutation = new GraphQLObjectType({
  name: "Mutation",
  fields:{
    addTask: {
      type: TaskType,
      args: {
        name: {type: new GraphQLNonNull(GraphQLString)},
        input: {type: new GraphQLNonNull(GraphQLList(GraphQLString))},
        command: {type: new GraphQLNonNull(GraphQLString)},
        output: {type: new GraphQLNonNull(GraphQLString)},
      },
      resolve(parentValue, args){
        return axios.post('http://localhost:3000/tasks', {
          name: args.name,
          input: args.input,
          command: args.command,
          output: args.output
        })
        .then(res => res.data);
      }
    }
  }
});

module.exports = new GraphQLSchema({
  query: RootQuery,
  mutation
});
