#include <iostream>
using namespace std;

enum Houses { Heart, Diamonds, Clubs, Spades };
enum Face {  Ace, Two, Three, Four, Five, Six, Seven, Eight, Nine, Ten, Jack, Queen, King };

// enum correspond to integers, hence we can use them as indices to access an array
const char *cardNames[] = {"A", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "J", "Q", "K" };
const char houseNames[] = {'H', 'D', 'C', 'S' };

struct Card
{
    Face m_face;
    Houses m_house;
};

ostream& operator<<(ostream& os, const Card& c )
{
    os << cardNames[c.m_face] << "-" << houseNames[c.m_house];
    return os;
};

struct Deck
{
    Card m_cards[52];

    Deck()
    {
        for (size_t h = 0; h < 4; h++)
            for (size_t f = 0; f < 13; f++) {
                m_cards[h*13 + f] = {(Face) f, (Houses) h};  // cast integer to enums
            }
    }

    void shuffle()
    {
        for (size_t i = 0; i < 500; ++i) {
            auto i1 = std::rand() % 52;
            auto i2 = std::rand() % 52;
            std::swap(m_cards[i1], m_cards[i2]);
        }

    }

};

ostream& operator<<(ostream& os, const Deck& deck)
{
    for (auto& c: deck.m_cards)
        os << c << ',';
    os << '\n';
    return os;
}


int main()
{
    Deck deck{};
    cout << deck << "\n";

    deck.shuffle();
    cout << deck << "\n";

    return 0;
}
