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

const BASE_URL = 'http://localhost:5000'

function getTaskByURL(relativeURL) {
  return axios.get(`${BASE_URL}${relativeURL}`)
    .then(res => res.data);
}

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
      resolve:(task) => task.connections.map(getTaskByURL)
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
        return axios.get('http://localhost:5000/tasks/'+args.id)
          .then(res => res.data);
      }
    },
    allTasks: {
      type: new GraphQLList(TaskType),
      resolve: () => getTaskByURL('/tasks/')
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
        return axios.post('http://localhost:5000/tasks', {
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
