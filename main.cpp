#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <memory>
#include <queue>

struct L1
{
    // Manhattan distance
    template <std::size_t Dimensions, typename Scalar>
    static Scalar manhattanDist(const std::array<Scalar, Dimensions> &location1,
                                const std::array<Scalar, Dimensions> &location2)
    {
        auto abs = [](Scalar v)
        { return v >= 0 ? v : -v; };
        Scalar dist = 0;
        for (std::size_t i = 0; i < Dimensions; i++)
        {
            dist += abs(location1[i] - location2[i]);
        }
        return dist;
    }
};

struct SquaredL2
{
    template <typename Scalar>
    static Scalar euclideanSqr(const std::vector<Scalar> &location1,
                               const std::vector<Scalar> &location2, std::size_t Dimensions)
    {
        auto sqr = [](Scalar v)
        { return v * v; };
        Scalar dist = 0;
        for (std::size_t i = 0; i < Dimensions; i++)
        {
            dist += sqr(location1[i] - location2[i]);
        }
        return dist;
    }
};

template <typename Scalar = double, std::size_t Dimensions = 2>
class KDTree
{
private:
    struct Node;
    std::vector<Node *> nodes; // Store pointers to nodes
    std::size_t bucketSize;

public:
    using Point = std::vector<Scalar>;

    KDTree() { bucketSize = 2 * (Dimensions); }

    void startKDTree(std::vector<Point> &points, std::size_t dimension)
    {
        nodes.reserve(points.size());
        buildKDTree(points, 0);
    }

    std::vector<Point> findKNearestNeighbors(const Point &queryPoint, std::size_t k)
    {
        std::priority_queue<std::pair<Scalar, Point>> knnQueue;

        // Distancias al cuadrado de los ejes
        std::vector<Scalar> a = {0, 0, 0};
        Scalar d = 0;

        // Comenzar la búsqueda KNN desde la raíz del árbol
        findKNN(nodes[0], queryPoint, k, a, d, knnQueue);

        // Extraer los k vecinos más cercanos de la cola de prioridad
        std::vector<Point> kNearestNeighbors;
        while (!knnQueue.empty())
        {
            kNearestNeighbors.push_back(knnQueue.top().second);
            knnQueue.pop();
        }

        // Los vecinos más cercanos se almacenan en orden inverso en la cola de prioridad, así que los invertimos aquí
        std::reverse(kNearestNeighbors.begin(), kNearestNeighbors.end());

        return kNearestNeighbors;
    }

    void findKNN(Node *node, const Point &queryPoint, std::size_t k, std::vector<Scalar> &a, Scalar &d, std::priority_queue<std::pair<Scalar, Point>> &knnQueue)
    {
        if (node->isLeaf)
        {
            for (const Point &point : node->bucket)
            {
                Scalar distance = SquaredL2::euclideanSqr(queryPoint, point, Dimensions);
                knnQueue.push(std::make_pair(distance, point));
            }
        }
        else
        {
            // Calcular la distancia entre el punto de consulta y el nodo actual
            Scalar distance = SquaredL2::euclideanSqr(queryPoint, node->location, Dimensions);

            // Realizar búsqueda en el hijo más cercano primero
            Node *nearestChild = (queryPoint[node->splitDimension] < node->location[node->splitDimension]) ? node->leftChild : node->rightChild;
            Node *farthestChild = (nearestChild == node->leftChild) ? node->rightChild : node->leftChild;

            findKNN(nearestChild, queryPoint, k, a, d, knnQueue); 

            Scalar sqrDistanceU = SquaredL2::euclideanSqr(queryPoint, node->location, node->splitDimension);
            d = d - a[node->splitDimension] + sqrDistanceU; 
            a[node->splitDimension] = sqrDistanceU;

            if (d < knnQueue.top().first)
            {
                findKNN(farthestChild, queryPoint, k, a, d, knnQueue);
            }
        }
    }

    std::vector<Point>
    findKNearestNeighborsWithBacktracking(const Point &queryPoint, std::size_t k)
    {
        std::vector<Point> primaryPath; // Almacena el camino primario
        Node *currentNode = nodes[0];   // Comenzar desde la raíz

        // Encontrar el camino primario y detenerse cuando lleguemos a la hoja que contiene el punto de consulta
        while (!currentNode->isLeaf)
        {
            primaryPath.push_back(currentNode->location);
            if (queryPoint[currentNode->splitDimension] < currentNode->location[currentNode->splitDimension])
            {
                currentNode = currentNode->leftChild;
            }
            else
            {
                currentNode = currentNode->rightChild;
            }
        }

        primaryPath.push_back(currentNode->location); // Agregar la hoja que contiene el punto de consulta

        std::priority_queue<std::pair<Scalar, Point>> knnQueue;

        // Realizar el algoritmo 1 (FindKNN) en la hoja que contiene el punto de consulta
        findKNN(currentNode, queryPoint, k, knnQueue);

        // Realizar el backtracking desde el camino primario
        for (int i = primaryPath.size() - 2; i >= 0; i--)
        {
            Node *node = nodes[i];
            Scalar axisDistance = std::abs(queryPoint[node->splitDimension] - node->location[node->splitDimension]);
            Scalar maxDistance = knnQueue.top().first;

            if (axisDistance < maxDistance)
            {
                // Si la distancia en la dimensión actual es menor que la máxima distancia actual, buscar en el otro hijo
                Node *otherChild = (queryPoint[node->splitDimension] < node->location[node->splitDimension]) ? node->rightChild : node->leftChild;
                findKNN(otherChild, queryPoint, k, knnQueue);
            }
        }

        // Extraer los k vecinos más cercanos de la cola de prioridad
        std::vector<Point> kNearestNeighbors;
        while (!knnQueue.empty())
        {
            kNearestNeighbors.push_back(knnQueue.top().second);
            knnQueue.pop();
        }

        // Los vecinos más cercanos se almacenan en orden inverso en la cola de prioridad, así que los invertimos aquí
        std::reverse(kNearestNeighbors.begin(), kNearestNeighbors.end());

        return kNearestNeighbors;
    }

    // Getters & setters

    Scalar getBucketSize() const { return bucketSize; }

    void setBucketSize(std::size_t bucketSize) { this->bucketSize = bucketSize; }

    ~KDTree()
    {
        for (Node *node : nodes)
        {
            delete node;
        }
    }

private:
    Node *buildKDTree(std::vector<Point> &points, std::size_t dimension)
    {
        if (points.size() <= bucketSize)
        {
            // Si encontramos un nodo hoja, guárdalo en el vector de nodos y usa max size_t para la dimensión de división
            std::size_t dim = std::numeric_limits<std::size_t>::max();
            return new Node(points, dim);
        }

        std::size_t lenDivided = points.size() / 2;

        // Ordenar puntos por la dimensión actual
        std::sort(points.begin(), points.end(), [dimension](const Point &a, const Point &b)
                  { return a[dimension] < b[dimension]; });

        Node *medianNode = new Node(points[lenDivided], dimension);

        std::vector<Point> leftPoints(points.begin(), points.begin() + lenDivided);
        std::vector<Point> rightPoints(points.begin() + lenDivided, points.end());

        nodes.push_back(medianNode);

        medianNode->leftChild = buildKDTree(leftPoints, (dimension + 1) % Dimensions);
        medianNode->rightChild = buildKDTree(rightPoints, (dimension + 1) % Dimensions);

        return medianNode;
    }

    struct Node
    {
        Point location;
        std::size_t splitDimension;
        Node *leftChild;
        Node *rightChild;
        std::vector<Point> bucket;
        bool isLeaf;

        Node(Point &location, std::size_t &splitDimension) : location(location), splitDimension(splitDimension), leftChild(nullptr), rightChild(nullptr), isLeaf(false) { bucket = {}; }
        Node(std::vector<Point> &locations, std::size_t &splitDimension) : location(Point()), splitDimension(splitDimension), leftChild(nullptr), rightChild(nullptr), bucket(locations), isLeaf(true) {}
    };
};

int main()
{
    std::vector<std::vector<double>> points =
        {
            {2, 4}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7, 2}, {1, 1}, {3, 3}, {6, 6}, {7, 8}, {2, 9}, {5, 2}, {8, 5}, {3, 7}, {4, 4}, {9, 3}, {1, 9}, {7, 7}, {5, 8}, {2, 7}, {4, 2}, {6, 4}, {9, 1}, {8, 9}, {2, 2}};

    std::vector<double> queryPoint = {2, 3};

    KDTree<double, 2> tree;

    tree.startKDTree(points, 0);

    // Realizar una búsqueda KNN utilizando el algoritmo 1
    std::size_t k = 10; // Número de vecinos más cercanos a buscar
    std::vector<std::vector<double>> result1 = tree.findKNearestNeighbors(queryPoint, k);

    // Realizar una búsqueda KNN utilizando el algoritmo 2 (con backtracking)
    // std::vector<std::vector<double>> result2 = tree.findKNearestNeighborsWithBacktracking(queryPoint, k);

    // Imprimir los resultados
    std::cout << "Resultados de la busqueda KNN (Algoritmo 1):\n";
    for (const auto &point : result1)
    {
        std::cout << "(" << point[0] << ", " << point[1] << ")\n";
    }

    /*  std::cout << "\nResultados de la búsqueda KNN (Algoritmo 2 - Backtracking):\n";
      for (const auto &point : result2)
      {
          std::cout << "(" << point[0] << ", " << point[1] << ")\n";
      }*/

    return 0;
}