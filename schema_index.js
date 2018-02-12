const fetch =require('node-fetch');
const axios = require('axios');
const {
  GraphQLID,
  GraphQLList,
  GraphQLNonNull,
  GraphQLObjectType,
  GraphQLSchema,
  GraphQLString,
} = require('graphql');
const {
  fromGlobalId,
  globalIdField,
  nodeDefinitions,
} = require('graphql-relay');

const {
  nodeField,
  nodeInterface,
} = nodeDefinitions(
  // A method that maps from a global id to an object
  (globalId, {loaders}) => {
    const {id, type} = fromGlobalId(globalId);
    if (type === 'Task') {
      return loaders.task.load(id);
    }
  },
  // A method that maps from an object to a type
  (obj) => {
    if (obj.hasOwnProperty('name')) {
      return TaskType;
    }
  }
);

const TaskType = new GraphQLObjectType({
  name: 'Task',
  description: 'A WDL Task',
  fields: () => ({
    id: globalIdField('Task'),
    name: {
      type: GraphQLString,
      description: 'task name',
      resolve: obj => obj.name,
    },
    input: {
      type: GraphQLList(GraphQLString),
      description: 'List of inputs',
      resolve: obj => obj.input,
    },
    command: {
      type: GraphQLString,
      description: 'The task command to be run',
      resolve: obj => obj.command,
    },
    output: {
      type: GraphQLString,
      description: 'Output of the task',
      resolve: obj => obj.output,
    },
    connections: {
      type: new GraphQLList(TaskType),
      description: 'Tasks that make a workflow',
      resolve: (obj, args, {loaders}) =>
        loaders.task.loadManyByURL(obj.connections),
    },
  }),
  interfaces: [nodeInterface],
});

const QueryType = new GraphQLObjectType({
  name: 'Query',
  description: 'The root of all... queries',
  fields: () => ({
    allTasks: {
      type: new GraphQLList(TaskType),
      description: 'All Tasks',
      resolve: (root, args, {loaders}) => loaders.task.loadAll(),
    },
    node: nodeField,
    task: {
      type: TaskType,
      args: {
        id: {type: new GraphQLNonNull(GraphQLID)},
      },
      resolve: (root, args, {loaders}) => loaders.task.load(args.id),
    },
  }),
});

module.export = new GraphQLSchema({
  query: QueryType,
});
