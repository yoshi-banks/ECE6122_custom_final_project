#include <arpa/inet.h>      // for inet_pton
#include <unistd.h>         // for close, read, write
#include <sys/socket.h>     // for socket, connect, send, recv
#include <netinet/in.h>     // for sockaddr_in
#include <cstring>          // for memset#include <iostream>
#include <iostream>
#include <thread>
#include <chrono>
#include "TerrainServerOpenTopo.hpp"

// Test client to verify server is working
void testClient(int port, double lat, double lon)
{
    std::this_thread::sleep_for(std::chrono::seconds(1)); // Wait for server to start

    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0)
    {
        std::cerr << "Client: Socket creation failed\n";
        return;
    }

    sockaddr_in serv_addr{};
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);

    if (inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr) <= 0)
    {
        std::cerr << "Client: Invalid address\n";
        close(sock);
        return;
    }

    if (connect(sock, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0)
    {
        std::cerr << "Client: Connection failed\n";
        close(sock);
        return;
    }

    std::string request = std::to_string(lat) + "," + std::to_string(lon);
    send(sock, request.c_str(), request.size(), 0);

    char buffer[1024] = {0};
    read(sock, buffer, sizeof(buffer));
    std::cout << "Client received: " << buffer;

    close(sock);
}

int main()
{
    std::cout << "=== TCP Server Test ===\n\n";

    TerrainServerOpenTopo server;
    int port = 8080;

    std::cout << "Starting server on port " << port << "...\n";
    std::cout << "Press Ctrl+C to stop the server.\n\n";

    // Launch a test client in a separate thread
    std::thread client_thread([port]() {
        std::cout << "Test client: Querying Grand Canyon South Rim (36.0544, -112.1401)...\n";
        testClient(port, 36.0544, -112.1401);

        std::this_thread::sleep_for(std::chrono::seconds(2));

        std::cout << "\nTest client: Querying Mount Everest (27.9881, 86.9250)...\n";
        testClient(port, 27.9881, 86.9250);
    });

    client_thread.detach();

    // Start the server (this will block)
    server.startServer(port);

    return 0;
}