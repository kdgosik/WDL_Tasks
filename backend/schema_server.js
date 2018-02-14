const axios = require('axios');
const {
  GraphQLObjectType,
  GraphQLString,
  GraphQLInt,
  GraphQLSchema,
  GraphQLList,
  GraphQLNonNull
} = require('graphql');


// This is where json-server is servinb
const BASE_URL = 'http://localhost:5000'

function getByURL(relativeURL) {
  return axios.get(`${BASE_URL}${relativeURL}`)
    .then(res => res.data);
}

// Defining Input Type
const InputType = new GraphQLObjectType({
  name:'Input',
  fields:() => ({
    id: {type: GraphQLInt},
    name: {type: GraphQLString},
    type: {type: GraphQLString}
  })
});


// Defining TaskType
const TaskType = new GraphQLObjectType({
  name:'Task',
  fields:() => ({
    id: {type: GraphQLInt},
    name: {type: GraphQLString},
    inputs: {
      type: GraphQLList(InputType),
      description: " Identifying inputs",
      resolve:(task) => task.inputs.map(getByURL)
    },
    command: {type: GraphQLString},
    output: {type: GraphQLString},
    connections: {
      type: GraphQLList(TaskType),
      description: "Finding related tasks",
      resolve:(task) => task.connections.map(getByURL)
    }
  })
});

// Root Query
const query = new GraphQLObjectType({
  name: 'Query',
  fields: {
    inputs: {
      type: InputType,
      args: {
        id:{type: GraphQLInt}
      },
      resolve(parentValue, args){
        return axios.get('http://localhost:5000/inputs/'+args.id)
          .then(res => res.data);
      }
    },
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
    allInputs: {
      type: new GraphQLList(TaskType),
      resolve: () => getTaskByURL('/inputs/')
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
  query,
  mutation
});
