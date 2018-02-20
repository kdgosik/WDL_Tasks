const axios = require('axios');
const {
  GraphQLObjectType,
  GraphQLString,
  GraphQLInt,
  GraphQLSchema,
  GraphQLList,
  GraphQLNonNull,
  GraphQLBoolean
} = require('graphql');


// This is where json-server is serving
const BASE_URL = 'http://localhost:5000'

function getByURL(relativeURL) {
  return axios.get(`${BASE_URL}${relativeURL}`)
    .then(res => res.data);
}

// Defining Input Type
const InputType = new GraphQLObjectType({
  name:'Input',
  fields:() => ({
    id: {type: GraphQLString},
    name: {type: GraphQLString},
    type: {type: GraphQLString}
  })
});


// Defining TaskType
const TaskType = new GraphQLObjectType({
  name:'Task',
  fields:() => ({
    id: {type: GraphQLString},
    name: {type: GraphQLString},
    inputs: {
      type: GraphQLList(InputType),
      description: " Identifying inputs",
      resolve:(obj) => obj.inputs.map(getByURL)
    },
    command: {type: GraphQLString},
    output: {type: GraphQLString},
    connections: {
      type: GraphQLList(TaskType),
      description: "Finding related tasks",
      resolve:(obj) => obj.connections.map(getByURL)
    }
  })
});


const ToolsType = new GraphQLObjectType({
  name: 'Tools',
  fields:() => ({
    id: {type: GraphQLString},
    organization: {type: GraphQLString},
    author: {type: GraphQLString},
    toolname: {type: GraphQLString},
    url: {type: GraphQLString},
    description: {type: GraphQLString},
    versions: {type: GraphQLList(GraphQLString)},
    verified_source: {type: GraphQLString},
    signed: {type: GraphQLString},
    meta_version: {
      type: GraphQLString,
      resolve: (obj, args, context) => obj['meta-version']
    },
    toolclass: {type: GraphQLList(GraphQLString)},
    contains: {type: GraphQLList(GraphQLString)},
    verified: {type: GraphQLBoolean}
  })
})

// Root Query
const query = new GraphQLObjectType({
  name: 'Query',
  fields: {
    inputs: {
      type: InputType,
      args: {
        id:{type: GraphQLString}
      },
      resolve(obj, args, context){
        return axios.get('http://localhost:5000/inputs/'+args.id)
          .then(res => res.data);
      }
    },
    task: {
      type: TaskType,
      args: {
        id:{type: GraphQLString}
      },
      resolve(obj, args, context){
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
    },
    allTools: {
      type: new GraphQLList(ToolsType),
      resolve: () => axios.get('https://api.firecloud.org/ga4gh/v1/tools').then(res => res.data)
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
      resolve(obj, args, context){
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
