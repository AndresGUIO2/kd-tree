#include <iostream>
#include <vector>
#include <algorithm>
#include <array>
#include <memory>
#include <queue>
#include <chrono>

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
        bucketSize = 2;
    };

    void startKDTree(std::vector<Point> &points, std::size_t dimension)
    {
        nodes.reserve(points.size());

        buildKDTree(points, 0);
    }

    void startKDTree(std::vector<Point> &Points)
    {
        nodes.reserve(Points.size());

        buildKDTree(Points);
    }

    std::vector<Point> findKNearestNeighbors(const Point &queryPoint, std::size_t k)
    {
        std::priority_queue<std::pair<Scalar, Point>> knnQueue;

        // Distancias al cuadrado de los ejes
        std::vector<Scalar> a(Dimensions, 0);
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
            return new Node(points, dim);
        }

        std::size_t medianIndex = points.size() / 2;

        // Ordenar puntos por la dimensión actual (introsort)
        std::sort(points.begin(), points.end(), [dimension](const Point &a, const Point &b)
                  { return a[dimension] < b[dimension]; });

        Node *medianNode = new Node(points[medianIndex], dimension);

        std::vector<Point> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<Point> rightPoints(points.begin() + medianIndex + 1, points.end());

        nodes.push_back(medianNode);

        medianNode->leftChild = buildKDTree(leftPoints, (dimension + 1) % Dimensions);
        medianNode->rightChild = buildKDTree(rightPoints, (dimension + 1) % Dimensions);

        return medianNode;
    }

    Node *buildKDTree(std::vector<Point> &points)
    {
        if (points.size() <= bucketSize)
        {
            std::size_t dim = std::numeric_limits<std::size_t>::max();
            return new Node(points, dim);
        }

        std::size_t dimension = selectSplitDimension(points);
        std::cout << "Dimension: " << dimension << std::endl;

        std::size_t medianIndex = points.size() / 2;

        std::sort(points.begin(), points.end(), [dimension](const Point &a, const Point &b)
                  { return a[dimension] < b[dimension]; });

        Node *medianNode = new Node(points[medianIndex], dimension);

        std::vector<Point> leftPoints(points.begin(), points.begin() + medianIndex);
        std::vector<Point> rightPoints(points.begin() + medianIndex + 1, points.end());

        nodes.push_back(medianNode);

        medianNode->leftChild = buildKDTree(leftPoints);
        medianNode->rightChild = buildKDTree(rightPoints);

        return medianNode;
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
    const int numVectores = 2000000;
    const int dimension = 4;

    // Crear un vector de vectores para almacenar los vectores de 4 dimensiones
    std::vector<std::vector<int>> vectores;

    // Generar 20000 vectores de 4 dimensiones con valores aleatorios
    for (int i = 0; i < numVectores; ++i)
    {
        std::vector<int> vector4D;
        for (int j = 0; j < dimension; ++j)
        {
            // Generar números aleatorios entre 1 y 100 para cada dimensión
            int valor = rand() % 1000000 + 1;
            vector4D.push_back(valor);
        }
        vectores.push_back(vector4D);
    }

    KDTree<int, dimension> kdTree;

    kdTree.startKDTree(vectores);

    std::vector<int> queryPoint = {58775, 205331, 31236, 47351};

    std::vector<std::vector<int>> kNearestNeighbors = kdTree.findKNearestNeighbors(queryPoint, 5);

    //mostrar los k vecinos

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < kNearestNeighbors.size(); ++i)
    {
        std::cout << "Vecino " << i + 1 << ": ";
        for (int j = 0; j < dimension; ++j)
        {
            std::cout << kNearestNeighbors[i][j] << " ";
        }
        std::cout << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Tiempo de ejecución: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " microsegundos" << std::endl;


    return 0;
}