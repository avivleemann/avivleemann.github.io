---
date: "2022-10-24"
sections:
- block: about.avatar
  content:
    text: null
    username: admin
    title: about me
  id: about
- block: features
  content:
    items:
    - icon: r-project
      icon_pack: fab
      name: R
    - icon: chart-line
      icon_pack: fas
      name: Statistics
    - icon: python
      icon_pack: fab
      name: python
    - icon: linux
      icon_pack: fab
      name: linux
    title: Skills
    
- block: portfolio
  content:
    # Choose how many pages you would like to display (0 = all pages)
    count: 0
    filters:
      folders:
        -projects
    title: Projects
    # Choose how many pages you would like to offset by
    offset: 0
    # Page order: descending (desc) or ascending (asc) date.
    order: desc
  
  design:
    columns: "1"
    flip_alt_rows: false
    # Choose a view for the listings:
    #   1 = List
    #   2 = Compact
    #   3 = Card
    #   4 = Citation (publication only)
    view: 3

type: landing

---


