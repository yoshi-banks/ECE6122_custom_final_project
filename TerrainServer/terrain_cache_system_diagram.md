```mermaid
flowchart TD
    Start([getTerrainData Called]) --> LoadCaches[Load All Cache Files on Startup]
    LoadCaches --> ExactMatch{Exact Cache<br/>Match?}
    
    ExactMatch -->|Yes| ReturnCache[Return Cached Data]
    ExactMatch -->|No| GenerateGrid[Generate Requested Grid Points]
    
    GenerateGrid --> CheckEachPoint[For Each Requested Point]
    CheckEachPoint --> WithinTolerance{Point Within<br/>Cached Bbox<br/>+ Tolerance?}
    
    WithinTolerance -->|Yes| Interpolate[Bilinear Interpolation<br/>from Cached Data]
    WithinTolerance -->|No| AddToMissing[Add to Missing List]
    
    Interpolate --> NextPoint{More<br/>Points?}
    AddToMissing --> NextPoint
    
    NextPoint -->|Yes| CheckEachPoint
    NextPoint -->|No| AnyMissing{Missing<br/>Points?}
    
    AnyMissing -->|Yes| FetchAPI[Fetch Missing Points<br/>from API in Batches]
    AnyMissing -->|No| Combine[Combine Interpolated Points]
    
    FetchAPI --> Combine
    Combine --> SaveNew[Save Complete Dataset<br/>to New Cache File]
    SaveNew --> AddToLoaded[Add to loadedCaches_]
    AddToLoaded --> ReturnData[Return Complete Data]
    
    ReturnCache --> End([End])
    ReturnData --> End
    
    style Start fill:#e1f5ff
    style End fill:#e1f5ff
    style ExactMatch fill:#fff4e1
    style WithinTolerance fill:#fff4e1
    style AnyMissing fill:#fff4e1
    style Interpolate fill:#d4edda
    style FetchAPI fill:#f8d7da
    style ReturnCache fill:#d4edda

