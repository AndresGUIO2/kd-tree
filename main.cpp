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
    std::vector<Node *> nodes;
    std::size_t bucketSize;

public:
    using Point = std::vector<Scalar>;

    KDTree()
    {
        bucketSize = nextPowerOf2(2 * (Dimensions));
        std::cout << "bucketSize: " << bucketSize << std::endl;
    };

    void startKDTree(std::vector<Point> &points, std::size_t dimension)
    {
        nodes.reserve(points.size());

        buildKDTree(points, 0);
    }

    void startAdaptiveKDTree(std::vector<Point> &Points, std::size_t dimension)
    {
        nodes.reserve(Points.size());

        Scalar splitDimension = selectSplitDimension(Points);

        std::cout << "Dimension de division: " << splitDimension << std::endl;

        buildAdaptiveKDTree(Points, splitDimension);
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
                if (knnQueue.size() == k)
                {
                    // Si la cola de prioridad está llena, comprueba si la distancia entre el punto de consulta y el punto actual es menor que la distancia del vecino más lejano
                    Scalar distance = SquaredL2::euclideanSqr(queryPoint, point, Dimensions);
                    if (distance < knnQueue.top().first)
                    {
                        // Si es así, elimina el vecino más lejano y añade el punto actual a la cola de prioridad
                        knnQueue.pop();
                        knnQueue.push(std::make_pair(distance, point));
                    }
                }
                else
                {
                    // Si la cola de prioridad no está llena, añade el punto actual a la cola de prioridad
                    Scalar distance = SquaredL2::euclideanSqr(queryPoint, point, Dimensions);
                    knnQueue.push(std::make_pair(distance, point));
                }
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

            if (d < knnQueue.top().first || knnQueue.size() < k)
            {
                findKNN(farthestChild, queryPoint, k, a, d, knnQueue);
            }

            if (distance < knnQueue.top().first)
            {
                // elimina el vecino más lejano y añade el punto actual a la cola de prioridad
                knnQueue.pop();
                knnQueue.push(std::make_pair(distance, node->location));
            }
        }
    }

    std::size_t nextPowerOf2(std::size_t n)
    {
        /*if (n && !(n & (n - 1)))
        {
            return n;
        }*/

        std::size_t count = 0;

        std::size_t original_n = n;

        while (n != 0)
        {
            n >>= 1;
            count += 1;
        }

        std::size_t result = 1 << count - 1;

        if (result < original_n)
        {
            return result << 1;
        }
        else
        {
            return result;
        }
    }

    int main()
    {
        std::cout << nextPowerOf2(1) << std::endl;  // Output: 1
        std::cout << nextPowerOf2(2) << std::endl;  // Output: 2
        std::cout << nextPowerOf2(3) << std::endl;  // Output: 4
        std::cout << nextPowerOf2(16) << std::endl; // Output: 16
        std::cout << nextPowerOf2(17) << std::endl; // Output: 32
        return 0;
    }

    // Getters & setters

    Scalar getBucketSize() const { return bucketSize; };

    void setBucketSize(std::size_t bucketSize) { this->bucketSize = bucketSize; };

    void setDimensions(std::size_t dimensions) { this->dimensions = dimensions; };

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
            std::cout << "Nodo hoja creado"
                      << "\n";
            return new Node(points, dim);
        }

        std::size_t medianIndex = points.size() / 2;

        //Ordenar puntos por la dimensión actual (introsort)
        std::sort(points.begin(), points.end(), [dimension](const Point &a, const Point &b)
                { return a[dimension] < b[dimension]; });

        Node *medianNode = new Node(points[medianIndex], dimension);

        std::cout << "Nodo raiz o interno: " << medianNode->location[0] << ", " << medianNode->location[1] << std::endl;

        std::vector<Point> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<Point> rightPoints(points.begin() + medianIndex + 1, points.end());

        nodes.push_back(medianNode);

        medianNode->leftChild = buildKDTree(leftPoints, (dimension + 1) % Dimensions);
        medianNode->rightChild = buildKDTree(rightPoints, (dimension + 1) % Dimensions);

        return medianNode;
    }

    Node *buildAdaptiveKDTree(std::vector<Point> &points, std::size_t dimension)
    {

        if (points.size() <= bucketSize)
        {
            std::size_t dim = std::numeric_limits<std::size_t>::max();
            std::cout << "Nodo hoja creado"
                      << "\n";
            return new Node(points, dim);
        }
    }

    int selectSplitDimension(const std::vector<Point> &points)
    {

        std::vector<Scalar> dispersion(Dimensions, 0.0);
        std::vector<Scalar> mean(Dimensions, 0.0);
        std::size_t n = points.size();

        // Calcular la suma de las coordenadas de cada dimensión
        for (const Point &point : points)
        {
            for (int dimension = 0; dimension < Dimensions; ++dimension)
            {
                mean[dimension] += point[dimension];
            }
        }
        // Calcular la media de las coordenadas de cada dimensión
        for (int dimension = 0; dimension < Dimensions; ++dimension)
        {
            mean[dimension] /= n;
        }

        // Calcular la "dispersión" (suma de diferencias absolutas) para cada dimensión
        for (const Point &point : points)
        {
            for (int dimension = 0; dimension < Dimensions; ++dimension)
            {
                dispersion[dimension] += (point[dimension] > mean[dimension] ? point[dimension] - mean[dimension] : mean[dimension] - point[dimension]);
            }
        }
        for (int dimension = 0; dimension < Dimensions; ++dimension)
        {
            dispersion[dimension] /= n;
        }

        // Encontrar la dimensión con la mayor dispersión
        int splitDimension = 0;
        Scalar maxDispersion = -std::numeric_limits<Scalar>::max();
        for (int dimension = 0; dimension < Dimensions; ++dimension)
        {
            if (dispersion[dimension] > maxDispersion)
            {
                maxDispersion = dispersion[dimension];
                splitDimension = dimension;
            }
        }

        std::cout << "Se ha seleccionado " << splitDimension << std::endl;

        return splitDimension;
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
            {2, 4}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7, 2}, {1, 1}, {8, 5}, {3, 7}, {4, 4}, {9, 3}, {1, 9}, {7, 7}, {5, 8}, {5, 7}};

    std::vector<double> queryPoint = {5, 4};

    KDTree<double, 2> tree;

    tree.startKDTree(points, 0);

    // Realizar una búsqueda KNN utilizando el algoritmo 1
    std::size_t k = 5; // Número de vecinos más cercanos a buscar
    std::vector<std::vector<double>> result1 = tree.findKNearestNeighbors(queryPoint, k);

    // Imprimir los resultados
    std::cout << "Resultados de la busqueda KNN (Algoritmo 1):\n";
    for (const auto &point : result1)
    {
        std::cout << "(" << point[0] << ", " << point[1] << ")\n";
    }

    return 0;
}