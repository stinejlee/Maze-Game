import java.util.ArrayList;
import java.util.Arrays;

import tester.*;
import javalib.impworld.*;
import java.awt.Color;
import javalib.worldimages.*;
import java.util.HashMap;
import java.util.Random;
import java.util.Deque;
import java.util.ArrayDeque;

// represents a vertex on a graph
class Vertex {
  Posn posn;
  ArrayList<Edge> outEdges;

  Vertex(Posn posn) {
    this.posn = posn;
    this.outEdges = new ArrayList<Edge>();
  }

  // adds the given edge to this Vertex's list of edges
  void addEdge(Edge e) {
    outEdges.add(e);
  }

  // checks if this Vertex is equal to the given object o
  // If object o is a Vertex, checks fields of that Vertex's posn to ensure they
  // are the equal
  @Override
  public boolean equals(Object o) {
    if (!(o instanceof Vertex)) {
      return false;
    }
    else {
      Vertex that = (Vertex) o;
      return this.posn.x == that.posn.x && this.posn.y == that.posn.y;
    }
  }

  @Override
  public int hashCode() {
    return this.posn.x * this.posn.x + this.posn.y * this.posn.y;
  }
}

// represents an edge between two vertices in a graph
class Edge {
  Vertex from;
  Vertex to;
  int weight;

  Edge(Vertex from, Vertex to, int weight) {
    this.from = from;
    this.to = to;
    this.weight = weight;
  }

}

// represents a graph made up of vertices and edges
class Graph {
  int width; // given in number of cells
  int height; // given in number of cells
  ArrayList<ArrayList<Vertex>> vertices = new ArrayList<ArrayList<Vertex>>();
  ArrayList<Edge> edges = new ArrayList<Edge>();
  Random random;

  Graph(int width, int height, int seed) {
    this.width = width;
    this.height = height;
    this.random = new Random(seed);
    this.makeVertices();
    this.makeAllEdges();
  }

  Graph(int width, int height) {
    this.width = width;
    this.height = height;
    this.random = new Random();
    this.makeVertices();
    this.makeAllEdges();
  }

  // this constructor is used solely for the purpose of testing some world methods
  Graph(int width, int height, int seed, ArrayList<ArrayList<Vertex>> vertices,
      ArrayList<Edge> edges) {
    this.width = width;
    this.height = height;
    this.random = new Random(seed);
    this.vertices = vertices;
    this.edges = edges;
  }

  // this constructor is used solely for the purpose of testing some world methods
  Graph(int width, int height, int seed, boolean foo) {
    this.width = width;
    this.height = height;
    this.random = new Random(seed);
  }

  // creates all the vertices for the graph
  void makeVertices() {
    for (int i = 0; i < this.width; i++) {
      ArrayList<Vertex> temp = new ArrayList<Vertex>();
      for (int n = 0; n < this.height; n++) {
        temp.add(new Vertex(new Posn(i, n)));
      }
      vertices.add(temp);
    }
  }

  // forms all of the edges between adjacent vertices
  void makeAllEdges() {
    for (int i = 0; i < this.width; i++) {
      for (int n = 0; n < this.height; n++) {
        if (i < this.width - 1) {
          edges.add(
              new Edge(vertices.get(i).get(n), vertices.get(i + 1).get(n), random.nextInt(1000)));
        }
        if (n < this.height - 1) {
          edges.add(
              new Edge(vertices.get(i).get(n), vertices.get(i).get(n + 1), random.nextInt(1000)));
        }
      }
    }
  }
}

// represents a maze
class Maze {
  int width;
  int height;
  Graph graph;
  HashMap<Vertex, Vertex> hmap = new HashMap<Vertex, Vertex>();
  ArrayList<Edge> edgesInTree = new ArrayList<Edge>();

  Maze(Graph graph) {
    this.width = graph.width;
    this.height = graph.height;
    this.graph = graph;
    this.startHash(graph.vertices);
    this.kruskal(graph.edges);
  }

  // this constructor is used solely for the purpose of testing some world methods
  Maze(Graph graph, boolean foo) {
    this.width = graph.width;
    this.height = graph.height;
    this.graph = graph;
  }

  // forms the hashmap of vertices
  void startHash(ArrayList<ArrayList<Vertex>> v) {
    for (int i = 0; i < width; i++) {
      for (int n = 0; n < height; n++) {
        hmap.put(v.get(i).get(n), v.get(i).get(n));
      }
    }
  }

  // sorts the list of edges then runs kruskals algorithm to create a tree
  void kruskal(ArrayList<Edge> e) {
    for (int i = 0; i < e.size(); i++) {
      for (int j = e.size() - 1; j > i; j--) {
        if (e.get(i).weight > e.get(j).weight) {
          Edge tmp = e.get(i);
          e.set(i, e.get(j));
          e.set(j, tmp);
        }
      }
    }

    int i = 0;
    while (i < e.size()) {
      if (!this.findRep(hmap, e.get(i).from).equals(this.findRep(hmap, e.get(i).to))) {
        edgesInTree.add(e.get(i));
        e.get(i).from.addEdge(e.get(i));
        e.get(i).to.addEdge(new Edge(e.get(i).to, e.get(i).from, 0));
        this.union(hmap, this.findRep(hmap, e.get(i).from), this.findRep(hmap, e.get(i).to));
      }
      i++;
    }
  }

  // finds the representative vertex of the given vertex based on a given hashmap
  Vertex findRep(HashMap<Vertex, Vertex> hmap, Vertex v) {
    return this.findRepHelp(hmap, v);
  }

  // if the linked vertex is the same as the given vertex, return this vertex
  // otherwise run this method again on the given vertex's linked vertex
  Vertex findRepHelp(HashMap<Vertex, Vertex> hmap, Vertex v) {
    if (hmap.get(v).equals(v)) {
      return v;
    }
    else {
      return this.findRep(hmap, hmap.get(v));
    }
  }

  // puts the given vertex 1 as the data for vertex 2
  void union(HashMap<Vertex, Vertex> hmap, Vertex v1, Vertex v2) {
    hmap.put(v2, v1);
  }

  // checks the next item in worklist if it is the end point of the graph and
  // if is calls reconstruct otherwise does nothing if the next item was already
  // seen
  // otherwise adds all neighbors of the next item in worklist to worklist
  boolean searchHelp(Vertex from, Vertex to, ICollection<Vertex> worklist,
      Deque<Vertex> alreadySeen, Deque<Vertex> finalPath) {
    // as long as the worklist isn't empty...
    // while (!worklist.isEmpty()) {
    Vertex next = worklist.remove();
    if (next.equals(to)) {
      alreadySeen.addFirst(next);
      finalPath.add(next);
      this.reconstruct(alreadySeen, next, finalPath);
      return true;
    }
    else if (alreadySeen.contains(next)) {
      // do nothing: we've already seen this one
    }
    else {
      // add all the neighbors of next to the worklist for further processing
      for (Edge e : next.outEdges) {
        worklist.add(e.to);
      }
      // add next to alreadySeen, since we're done with it
      alreadySeen.addFirst(next);
    }
    // }
    // we haven't found the vertex, and there are no more to try
    return false;
  }

  // constructs the solution path through the maze
  void reconstruct(Deque<Vertex> alreadySeen, Vertex current, Deque<Vertex> finalPath) {
    for (Vertex v : alreadySeen) {
      for (Edge e : current.outEdges) {
        if (e.to.equals(v)) {
          current = e.to;
          finalPath.add(v);
        }
      }
    }
  }

}

// represents a world that is a maze game
class MazeGame extends World {

  Maze maze;
  int width;
  int height;
  int cellSize;
  boolean dfs = false;
  boolean bfs = false;
  boolean recon = false;
  ICollection<Vertex> worklist;
  Deque<Vertex> alreadySeen = new ArrayDeque<Vertex>();
  Deque<Vertex> finalPath = new ArrayDeque<Vertex>();
  Deque<Vertex> drawPath = new ArrayDeque<Vertex>();
  Vertex start;
  Vertex end;

  MazeGame(Maze maze) {
    this.maze = maze;
    if (maze.width >= maze.height) {
      this.width = 600;
      this.height = (int) (600 * ((double) maze.height / (double) maze.width));
      this.cellSize = 600 / maze.width;
    }
    else {
      this.width = (int) (600 * ((double) maze.width / (double) maze.height));
      this.height = 600;
      this.cellSize = 600 / maze.height;
    }
    this.start = this.maze.graph.vertices.get(0).get(0);
    this.end = this.maze.graph.vertices.get(this.maze.graph.width - 1)
        .get(this.maze.graph.height - 1);
  }

  // this constructor takes the desired cell dimensions
  MazeGame(int x, int y) {
    Graph graph = new Graph(x, y);
    Maze maze = new Maze(graph);
    this.maze = maze;
    if (maze.width >= maze.height) {
      this.width = 600;
      this.height = (int) (600 * ((double) maze.height / (double) maze.width));
      this.cellSize = 600 / maze.width;
    }
    else {
      this.width = (int) (600 * ((double) maze.width / (double) maze.height));
      this.height = 600;
      this.cellSize = 600 / maze.height;
    }
    this.start = this.maze.graph.vertices.get(0).get(0);
    this.end = this.maze.graph.vertices.get(this.maze.graph.width - 1)
        .get(this.maze.graph.height - 1);
  }

  // handles creating the image of the maze
  @Override
  public WorldScene makeScene() {
    WorldScene scene = new WorldScene(this.width, this.height);
    scene.placeImageXY(
        new RectangleImage(this.width, this.height, OutlineMode.SOLID, new Color(217, 189, 255)),
        this.width / 2, this.height / 2);
    scene.placeImageXY(new RectangleImage(this.cellSize, this.cellSize, OutlineMode.SOLID,
        new Color(208, 239, 255)), this.cellSize / 2, this.cellSize / 2);
    scene.placeImageXY(
        new RectangleImage(this.cellSize, this.cellSize, OutlineMode.SOLID,
            new Color(255, 184, 173)),
        this.cellSize * this.maze.width - this.cellSize / 2,
        this.cellSize * this.maze.height - this.cellSize / 2);
    this.drawWalls(scene, this.maze.graph.vertices, this.maze.edgesInTree);
    for (Vertex v : alreadySeen) {
      scene.placeImageXY(
          new RectangleImage((int) (this.cellSize * 0.8), (int) (this.cellSize * 0.8),
              OutlineMode.SOLID, Color.lightGray),
          v.posn.x * cellSize + cellSize / 2, v.posn.y * cellSize + cellSize / 2);
    }
    for (Vertex v : drawPath) {
      scene.placeImageXY(
          new RectangleImage((int) (this.cellSize * 0.8), (int) (this.cellSize * 0.8),
              OutlineMode.SOLID, Color.white),
          v.posn.x * cellSize + cellSize / 2, v.posn.y * cellSize + cellSize / 2);
    }

    return scene;
  }

  // handles the animation of the maze
  @Override
  public void onTick() {
    if (this.bfs) {
      this.bfs = !this.maze.searchHelp(start, end, worklist, alreadySeen, finalPath);
      this.recon = !this.bfs;
    }
    else if (this.dfs) {
      this.dfs = !this.maze.searchHelp(start, end, worklist, alreadySeen, finalPath);
      this.recon = !this.dfs;
    }
    if (this.recon) {
      if (finalPath.size() > 0) {
        Vertex v = finalPath.remove();
        drawPath.add(v);
      }
      else {
        this.recon = false;
      }
    }
  }

  // handles keys b and d to perform breadth first searches and depth first
  // searches respectively
  @Override
  public void onKeyEvent(String key) {
    if (key.equals("b")) {
      this.bfs = true;
      this.dfs = false;
      worklist = new Queue<Vertex>();
      worklist.add(this.start);
    }
    else if (key.equals("d")) {
      this.dfs = true;
      this.bfs = false;
      worklist = new Stack<Vertex>();
      worklist.add(this.start);
    }
  }

  // draws the walls of the maze
  void drawWalls(WorldScene scene, ArrayList<ArrayList<Vertex>> v, ArrayList<Edge> e) {
    for (int i = 0; i < v.size(); i++) {
      for (int n = 0; n < v.get(0).size(); n++) {
        if (i < this.width - 1) { // right wall
          scene.placeImageXY(new LineImage(new Posn(0, this.cellSize), Color.black),
              (i + 1) * this.cellSize, n * this.cellSize + this.cellSize / 2);
        }
        if (n < this.height - 1) { // bottom wall
          scene.placeImageXY(new LineImage(new Posn(this.cellSize, 0), Color.black),
              i * this.cellSize + this.cellSize / 2, (n + 1) * this.cellSize);
        }
      }
    }

    for (int i = 0; i < e.size(); i++) {
      if (e.get(i).from.posn.x == e.get(i).to.posn.x) {
        scene.placeImageXY(new LineImage(new Posn(this.cellSize, 0), new Color(217, 189, 255)),
            e.get(i).from.posn.x * this.cellSize + this.cellSize / 2,
            (e.get(i).from.posn.y * this.cellSize + e.get(i).to.posn.y * this.cellSize) / 2
                + this.cellSize / 2);
      }
      else if (e.get(i).from.posn.y == e.get(i).to.posn.y) {
        scene.placeImageXY(new LineImage(new Posn(0, this.cellSize), new Color(217, 189, 255)),
            (e.get(i).from.posn.x * this.cellSize + e.get(i).to.posn.x * this.cellSize) / 2
                + this.cellSize / 2,
            e.get(i).from.posn.y * this.cellSize + this.cellSize / 2);
      }
    }
  }
}

interface ICollection<T> {
  // is this collection empty?
  boolean isEmpty();

  // EFFECT: adds the item to the collection
  void add(T item);

  // returns the first item in the collection
  // EFFECT: removes the first item
  T remove();
}

class Stack<T> implements ICollection<T> {
  Deque<T> contents;

  Stack() {
    this.contents = new ArrayDeque<T>();
  }

  // is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  // EFFECT: adds the given item to the stack
  public void add(T item) {
    this.contents.addFirst(item);
  }

  // returns the first item in the collection
  // EFFECT: removes the first item
  public T remove() {
    return this.contents.removeFirst();
  }
}

class Queue<T> implements ICollection<T> {
  Deque<T> contents;

  Queue() {
    this.contents = new ArrayDeque<T>();
  }

  // is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  // EFFECT: adds the given item to the stack
  public void add(T item) {
    this.contents.addLast(item);
  }

  // returns the first item in the collection
  // EFFECT: removes the first item
  public T remove() {
    return this.contents.removeFirst();
  }
}

class ExamplesMazeWorld {
  ExamplesMazeWorld() {
  }

  MazeGame game0;

  Graph testGraph;
  Maze testMaze;
  MazeGame game;

  Graph testGraph2;
  Maze testMaze2;
  MazeGame game2;

  WorldScene scene;

  Graph testGraph3;
  Maze testMaze3;
  MazeGame game3;

  MazeGame game4;

  WorldScene scene2;

  WorldScene sceneWalls;

  Vertex a;
  Vertex a1;
  Vertex b;
  Vertex c;
  Vertex d;

  Edge ab;
  Edge ac;
  Edge bd;
  Edge cd;

  ArrayList<Edge> allEdges;
  ArrayList<Edge> allEdgesWrong;
  ArrayList<Edge> kruskalEdges;
  ArrayList<ArrayList<Vertex>> allVertices;
  Graph abcdGraph;
  Graph abcdGraph2;
  Graph mtGraph;
  Maze abcdMaze;
  Maze abcdMaze2;
  Maze mtMaze;

  HashMap<Vertex, Vertex> hash;
  HashMap<Vertex, Vertex> kruskal;
  HashMap<Vertex, Vertex> union;

  ICollection<Vertex> worklist1;
  ICollection<Vertex> worklist2;
  ICollection<Vertex> worklist3;
  ICollection<Vertex> worklist4;
  ICollection<Vertex> worklist5;

  Deque<Vertex> alreadySeen;
  Deque<Vertex> finalPath;
  Deque<Vertex> abcdFinalPath;

  void initData() {
    game0 = new MazeGame(100, 60);

    testGraph = new Graph(15, 15, 2);
    testMaze = new Maze(testGraph);
    game = new MazeGame(testMaze);

    testGraph2 = new Graph(3, 3, -1);
    testMaze2 = new Maze(testGraph2);
    game2 = new MazeGame(testMaze2);
    scene = new WorldScene(600, 600);

    testGraph3 = new Graph(6, 3, -1);
    testMaze3 = new Maze(testGraph3);
    game3 = new MazeGame(testMaze3);
    scene2 = new WorldScene(600, 300);

    a = new Vertex(new Posn(0, 0));
    a1 = new Vertex(new Posn(0, 0));
    b = new Vertex(new Posn(0, 1));
    c = new Vertex(new Posn(1, 0));
    d = new Vertex(new Posn(1, 1));

    ab = new Edge(a, b, 225);
    ac = new Edge(a, c, 913);
    bd = new Edge(b, d, 579);
    cd = new Edge(c, d, 439);

    allEdges = new ArrayList<Edge>(Arrays.asList(ac, ab, bd, cd));
    allEdgesWrong = new ArrayList<Edge>(Arrays.asList(cd, bd, ac, ab));
    kruskalEdges = new ArrayList<Edge>(Arrays.asList(ab, cd, bd));

    allVertices = new ArrayList<ArrayList<Vertex>>(Arrays.asList(
        new ArrayList<Vertex>(Arrays.asList(a, b)), new ArrayList<Vertex>(Arrays.asList(c, d))));

    abcdGraph = new Graph(2, 2, -1, allVertices, allEdges);
    abcdGraph2 = new Graph(2, 2, -1);
    mtGraph = new Graph(2, 2, -1, false);

    abcdMaze = new Maze(abcdGraph, false);
    abcdMaze2 = new Maze(abcdGraph2);
    mtMaze = new Maze(mtGraph, false);

    game4 = new MazeGame(abcdMaze2);

    hash = new HashMap<Vertex, Vertex>();
    hash.put(a, a);
    hash.put(b, b);
    hash.put(c, c);
    hash.put(d, d);
    kruskal = new HashMap<Vertex, Vertex>();
    kruskal.put(a, a);
    kruskal.put(b, b);
    kruskal.put(c, c);
    kruskal.put(d, d);
    kruskal.put(b, a);
    kruskal.put(d, c);
    kruskal.put(c, a);
    union = new HashMap<Vertex, Vertex>();
    union.put(a, a);
    union.put(b, b);
    union.put(c, c);
    union.put(d, d);
    union.put(b, a);

    sceneWalls = new WorldScene(600, 600);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (0 + 1) * 200,
        0 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (0 + 1) * 200,
        1 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (0 + 1) * 200,
        2 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (1 + 1) * 200,
        0 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (1 + 1) * 200,
        1 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), Color.black), (1 + 1) * 200,
        2 * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 0 * 200 + 200 / 2,
        (0 + 1) * 200);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 1 * 200 + 200 / 2,
        (0 + 1) * 200);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 2 * 200 + 200 / 2,
        (0 + 1) * 200);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 0 * 200 + 200 / 2,
        (1 + 1) * 200);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 1 * 200 + 200 / 2,
        (1 + 1) * 200);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), Color.black), 2 * 200 + 200 / 2,
        (1 + 1) * 200);

    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), new Color(217, 189, 255)),
        game2.maze.edgesInTree.get(0).from.posn.x * 200 + 200 / 2,
        (game2.maze.edgesInTree.get(0).from.posn.y * 200
            + game2.maze.edgesInTree.get(0).to.posn.y * 200) / 2 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), new Color(217, 189, 255)),
        (game2.maze.edgesInTree.get(1).from.posn.x * 200
            + game2.maze.edgesInTree.get(1).to.posn.x * 200) / 2 + 200 / 2,
        game2.maze.edgesInTree.get(1).from.posn.y * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), new Color(217, 189, 255)),
        game2.maze.edgesInTree.get(2).from.posn.x * 200 + 200 / 2,
        (game2.maze.edgesInTree.get(2).from.posn.y * 200
            + game2.maze.edgesInTree.get(2).to.posn.y * 200) / 2 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), new Color(217, 189, 255)),
        game2.maze.edgesInTree.get(3).from.posn.x * 200 + 200 / 2,
        (game2.maze.edgesInTree.get(3).from.posn.y * 200
            + game2.maze.edgesInTree.get(3).to.posn.y * 200) / 2 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), new Color(217, 189, 255)),
        (game2.maze.edgesInTree.get(4).from.posn.x * 200
            + game2.maze.edgesInTree.get(4).to.posn.x * 200) / 2 + 200 / 2,
        game2.maze.edgesInTree.get(4).from.posn.y * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(200, 0), new Color(217, 189, 255)),
        game2.maze.edgesInTree.get(5).from.posn.x * 200 + 200 / 2,
        (game2.maze.edgesInTree.get(5).from.posn.y * 200
            + game2.maze.edgesInTree.get(5).to.posn.y * 200) / 2 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), new Color(217, 189, 255)),
        (game2.maze.edgesInTree.get(6).from.posn.x * 200
            + game2.maze.edgesInTree.get(6).to.posn.x * 200) / 2 + 200 / 2,
        game2.maze.edgesInTree.get(6).from.posn.y * 200 + 200 / 2);
    sceneWalls.placeImageXY(new LineImage(new Posn(0, 200), new Color(217, 189, 255)),
        (game2.maze.edgesInTree.get(7).from.posn.x * 200
            + game2.maze.edgesInTree.get(7).to.posn.x * 200) / 2 + 200 / 2,
        game2.maze.edgesInTree.get(7).from.posn.y * 200 + 200 / 2);

    worklist1 = new Queue<Vertex>();
    worklist1.add(a);
    worklist2 = new Stack<Vertex>();
    worklist2.add(a);
    worklist3 = new Stack<Vertex>();
    worklist3.add(d);
    worklist4 = new Queue<Vertex>();
    worklist5 = new Stack<Vertex>();
    alreadySeen = new ArrayDeque<Vertex>();
    alreadySeen.add(a);
    alreadySeen.add(b);
    alreadySeen.add(c);
    alreadySeen.add(d);
    finalPath = new ArrayDeque<Vertex>();
    abcdFinalPath = new ArrayDeque<Vertex>();
    abcdFinalPath.add(b);
    abcdFinalPath.add(d);

  }

  // test method for the maze game world
  void testMazeGameWorld(Tester t) {
    this.initData();

    game0.bigBang(600, 600, 0.01);
  }

  // tests the addEdge method
  void testAddEdge(Tester t) {
    this.initData();

    t.checkExpect(a.outEdges, new ArrayList<Edge>());
    a.addEdge(ab);
    t.checkExpect(a.outEdges, new ArrayList<Edge>(Arrays.asList(ab)));
    a.addEdge(ac);
    t.checkExpect(a.outEdges, new ArrayList<Edge>(Arrays.asList(ab, ac)));
  }

  // tests the equals method
  void testEquals(Tester t) {
    this.initData();

    t.checkExpect(this.a.equals(b), false);
    t.checkExpect(this.a.equals(a1), true);
  }

  void testHashCode(Tester t) {
    this.initData();
    t.checkExpect(this.a.hashCode(), 0);
    t.checkExpect(this.b.hashCode(), 1);
    t.checkExpect(this.c.hashCode(), 1);
    t.checkExpect(this.d.hashCode(), 2);
  }

  // tests the makeVertices method
  void testMakeVertices(Tester t) {
    this.initData();

    t.checkExpect(this.mtGraph.vertices, new ArrayList<ArrayList<Vertex>>());
    this.mtGraph.makeVertices();
    t.checkExpect(this.mtGraph.vertices, this.allVertices);

  }

  // tests the makeAllEdges method
  void testMakeAllEdges(Tester t) {
    this.initData();

    this.mtGraph.makeVertices();
    t.checkExpect(this.mtGraph.edges, new ArrayList<Edge>());
    this.mtGraph.makeAllEdges();
    t.checkExpect(this.mtGraph.edges, this.allEdges);
  }

  // tests the startHash method
  void testStartHash(Tester t) {
    this.initData();

    t.checkExpect(this.abcdMaze.hmap, new HashMap<Vertex, Vertex>());
    this.abcdMaze.startHash(allVertices);
    t.checkExpect(this.abcdMaze.hmap, this.hash);

  }

  // tests the kruskal method
  void testKruskal(Tester t) {
    this.initData();

    this.abcdMaze.startHash(this.allVertices);
    t.checkExpect(this.allEdgesWrong.get(0).weight, 439);
    t.checkExpect(this.allEdgesWrong.get(1).weight, 579);
    t.checkExpect(this.abcdMaze.edgesInTree, new ArrayList<Edge>());
    t.checkExpect(this.abcdMaze.hmap, this.hash);
    this.abcdMaze.kruskal(this.allEdgesWrong);
    t.checkExpect(this.allEdgesWrong.get(0).weight, 225);
    t.checkExpect(this.allEdgesWrong.get(1).weight, 439);
    t.checkExpect(this.abcdMaze.edgesInTree, this.kruskalEdges);
    t.checkExpect(this.abcdMaze.hmap, this.kruskal);
  }

  // tests the findRep method
  void testFindRep(Tester t) {
    this.initData();

    t.checkExpect(this.abcdMaze.findRep(this.kruskal, this.a), this.a);
    t.checkExpect(this.abcdMaze.findRep(this.kruskal, this.d), this.a);

  }

  // tests the findRepHelp method
  void testFindRepHelp(Tester t) {
    this.initData();

    t.checkExpect(this.abcdMaze.findRepHelp(this.kruskal, this.a), this.a);
    t.checkExpect(this.abcdMaze.findRepHelp(this.kruskal, this.d), this.a);
  }

  // tests the union method
  void testUnion(Tester t) {
    this.initData();

    this.abcdMaze.startHash(this.allVertices);
    t.checkExpect(this.abcdMaze.hmap, this.hash);
    this.abcdMaze.union(this.abcdMaze.hmap, a, b);
    t.checkExpect(this.abcdMaze.hmap, this.union);
    this.abcdMaze.union(this.abcdMaze.hmap, c, d);
    this.abcdMaze.union(this.abcdMaze.hmap, a, c);
    t.checkExpect(this.abcdMaze.hmap, this.kruskal);

  }

  // tests the searchHelp method
  boolean testSearchHelp(Tester t) {
    this.initData();

    return t.checkExpect(this.abcdMaze.searchHelp(a, d, worklist1, alreadySeen, finalPath), false)
        && t.checkExpect(this.abcdMaze.searchHelp(a, d, worklist2, alreadySeen, finalPath), false)
        && t.checkExpect(this.abcdMaze.searchHelp(a, d, worklist3, alreadySeen, finalPath), true);
  }

  // tests the reconstruct method
  void testReconstruct(Tester t) {
    this.initData();

    t.checkExpect(finalPath, new ArrayDeque<Vertex>());
    this.abcdMaze2.reconstruct(alreadySeen, this.abcdMaze2.graph.vertices
        .get(this.abcdMaze2.graph.width - 1).get(this.abcdMaze2.graph.height - 1), finalPath);
    t.checkExpect(finalPath, abcdFinalPath);

  }

  // tests the makeScene method
  void testMakeScene(Tester t) {
    this.initData();

    scene.placeImageXY(new RectangleImage(600, 600, OutlineMode.SOLID, new Color(217, 189, 255)),
        300, 300);
    scene.placeImageXY(new RectangleImage(200, 200, OutlineMode.SOLID, new Color(208, 239, 255)),
        200 / 2, 200 / 2);
    scene.placeImageXY(new RectangleImage(200, 200, OutlineMode.SOLID, new Color(255, 184, 173)),
        200 * 3 - 200 / 2, 200 * 3 - 200 / 2);
    this.game2.drawWalls(this.scene, this.game2.maze.graph.vertices, this.game2.maze.edgesInTree);

    t.checkExpect(this.game2.makeScene(), this.scene);

    scene2.placeImageXY(new RectangleImage(600, 300, OutlineMode.SOLID, new Color(217, 189, 255)),
        600 / 2, 300 / 2);
    scene2.placeImageXY(new RectangleImage(100, 100, OutlineMode.SOLID, new Color(208, 239, 255)),
        100 / 2, 100 / 2);
    scene2.placeImageXY(new RectangleImage(100, 100, OutlineMode.SOLID, new Color(255, 184, 173)),
        100 * 6 - 100 / 2, 100 * 3 - 100 / 2);
    this.game3.drawWalls(this.scene2, this.game3.maze.graph.vertices, this.game3.maze.edgesInTree);

    t.checkExpect(this.game3.makeScene(), this.scene2);

  }

  // test method for onTick
  void testOnTick(Tester t) {
    this.initData();

    t.checkExpect(this.game4.bfs, false);
    t.checkExpect(this.game4.recon, false);
    this.game4.onKeyEvent("b");
    t.checkExpect(this.game4.bfs, true);
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.bfs, true);
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.bfs, false);
    t.checkExpect(this.game4.recon, true);
    this.game4.onTick();
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.recon, false);

    this.initData();

    t.checkExpect(this.game4.dfs, false);
    t.checkExpect(this.game4.recon, false);
    this.game4.onKeyEvent("d");
    t.checkExpect(this.game4.dfs, true);
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.dfs, true);
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.dfs, false);
    t.checkExpect(this.game4.recon, true);
    this.game4.onTick();
    this.game4.onTick();
    this.game4.onTick();
    t.checkExpect(this.game4.recon, false);
  }

  // test method for onKeyEvent
  void testOnKeyEvent(Tester t) {
    this.initData();

    t.checkExpect(this.game.bfs, false);
    t.checkExpect(this.game.dfs, false);
    this.game.onKeyEvent("b");
    t.checkExpect(this.game.bfs, true);
    t.checkExpect(this.game.dfs, false);
    this.game.onKeyEvent("d");
    t.checkExpect(this.game.bfs, false);
    t.checkExpect(this.game.dfs, true);
  }

  // test method for the draw walls method
  void testDrawWalls(Tester t) {
    this.initData();

    t.checkExpect(this.scene, new WorldScene(600, 600));
    this.game2.drawWalls(this.scene, this.game2.maze.graph.vertices, this.game2.maze.edgesInTree);
    t.checkExpect(this.scene, this.sceneWalls);

  }

  // test method for isEmpty
  void testIsEmpty(Tester t) {
    this.initData();

    t.checkExpect(worklist1.isEmpty(), false);
    t.checkExpect(worklist2.isEmpty(), false);
    t.checkExpect(worklist4.isEmpty(), true);
    t.checkExpect(worklist5.isEmpty(), true);
  }

  // test method for add
  void testAdd(Tester t) {
    this.initData();

    t.checkExpect(worklist4.isEmpty(), true);
    t.checkExpect(worklist5.isEmpty(), true);
    worklist4.add(a);
    worklist5.add(a);
    t.checkExpect(worklist4, worklist1);
    t.checkExpect(worklist5, worklist2);

  }

  // test method for remove
  void testRemove(Tester t) {
    this.initData();

    t.checkExpect(worklist1.remove(), a);
    t.checkExpect(worklist2.remove(), a);

  }
}