import React from 'react'
import { gql, graphql } from 'react-apollo'
//import { Link } from 'react-router-dom'
const Link = require('react-router-dom');

const query = gql`{
  allTasks {
    id
    name
  }
}`

class HomeView extends React.Component {
  render() {
    let { data } = this.props
    if (data.loading) {
      return <div>Loading...</div>
    }
    return (
      <div>
        {data.allTasks.map((item, index) => (
          <p key={item.id}>
            <Link to={`/tasks/${item.id}/`}>
              {item.name}
            </Link>
          </p>
        ))}
      </div>
    )
  }
}

HomeView = graphql(query)(HomeView)
// I dont think this is working.  Keeps saing property is undefined above when trying to map
export default HomeView
