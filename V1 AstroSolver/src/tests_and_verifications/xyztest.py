def delete_linked_list(head, node):
    """
    Delete a node from a singly linked list.
    """

    if head is None:
        return None

    if head == node:
        return head.next

    prev = head
    curr = head.next
    while curr is not None:
        if curr == node:
            prev.next = curr.next
            return head
        prev = curr
        curr = curr.next

    return head