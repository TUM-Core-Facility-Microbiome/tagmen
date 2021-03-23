# This code is contributed by Atul Kumar (www.facebook.com/atul.kr.007)
# From https://www.geeksforgeeks.org/trie-insert-and-search/
# Adapted for DNA alphabets

class TrieNode(object):
    # Trie node class
    def __init__(self):
        self.children = [None] * 4

        # values is set for leaves
        self.values = []

    def is_leaf(self):
        return len(self.values) != 0


class Trie(object):
    # Trie data structure class
    def __init__(self):
        self.root = self._get_new_node()

    @staticmethod
    def _get_new_node():
        # Returns new trie node (initialized to NULLs)
        return TrieNode()

    @staticmethod
    def _char_to_index(ch):
        indices = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
        # private helper function
        # Converts key current character into index
        # use only 'a' through 'z' and lower case
        return indices.get(ch, -1)

    def insert(self, key, value):
        # If not present, inserts key into trie
        # If the key is prefix of trie node,
        # just marks leaf node
        p_crawl = self.root
        length = len(key)
        for level in range(length):
            index = self._char_to_index(key[level])

            # if current character is not present
            if not p_crawl.children[index]:
                p_crawl.children[index] = self._get_new_node()
            p_crawl = p_crawl.children[index]

        # mark last node as leaf
        p_crawl.values.append(value)
        return p_crawl

    def search(self, key):
        # Search key in the trie
        # Returns true if key presents
        # in trie, else false
        p_crawl = self.root
        length = len(key)
        for level in range(length):
            index = self._char_to_index(key[level])
            if not p_crawl.children[index]:
                return None
            p_crawl = p_crawl.children[index]

        if p_crawl is not None and p_crawl.is_leaf():
            return p_crawl
        else:
            return None

    def search_longest_matching_prefix(self, key):
        p_crawl = self.root
        length = len(key)
        for level in range(length):
            if p_crawl.is_leaf():
                return True, key[:level], p_crawl
            index = self._char_to_index(key[level])
            if not p_crawl.children[index]:
                return False, key[:level], p_crawl
            p_crawl = p_crawl.children[index]

        return True, key, p_crawl
