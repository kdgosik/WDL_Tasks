const DataLoader = require('dataloader');

const express = require('express');
const fetch = require('node-fetch');
const graphqlHTTP = require('express-graphql');
const schema = require('./schema_index.js');

const BASE_URL = 'http://localhost:3000';

function getJSONFromRelativeURL(relativeURL) {
  return fetch(`${BASE_URL}${relativeURL}`)
    .then(res => res.json());
}

function getTasks() {
  return getJSONFromRelativeURL('/tasks/')
    .then(json => json.tasks);
}

function getTask(id) {
  return getTaskByURL(`/tasks/${id}/`);
}

function getTasksByURL(relativeURL) {
  return getJSONFromRelativeURL(relativeURL)
    .then(json => json.tasks);
}

const app = express();

app.use(graphqlHTTP(req => {
  const cacheMap = new Map();
  const tasksLoader =
    new DataLoader(keys => Promise.all(keys.map(getTasks)), {cacheMap});
  const taskLoader =
    new DataLoader(keys => Promise.all(keys.map(getTask)), {
      cacheKeyFn: key => `/tasks/${key}/`,
      cacheMap,
    });
  const taskByURLLoader =
    new DataLoader(keys => Promise.all(keys.map(getTaskByURL)), {cacheMap});
  taskLoader.loadAll = tasksLoader.load.bind(taskLoader, '__all__');
  taskLoader.loadByURL = taskByURLLoader.load.bind(taskByURLLoader);
  taskLoader.loadManyByURL =
    taskByURLLoader.loadMany.bind(taskByURLLoader);
  const loaders = {task: taskLoader};
  return {
    context: {loaders},
    graphiql: true,
    schema,
  };
}));

app.listen(
  4000,
  () => console.log('GraphQL Server running at http://localhost:4000')
);
