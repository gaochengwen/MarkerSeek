function bySpecies(items, speciesGetter) {
  return items.reduce((groups, item) => {
    const species = speciesGetter(item) || "Unknown";
    if (!groups[species]) {
      groups[species] = [];
    }
    groups[species].push(item);
    return groups;
  }, {});
}

const PLOT_MARGIN = {l: 50, r: 30, t: 30, b: 50};
const PLOT_FONT = {family: "Inter, system-ui, sans-serif", size: 12};
const GRID_COLOR = "#e6e6e6";
const LEGEND = {
  x: 1,
  y: 1,
  xanchor: "right",
  yanchor: "top",
  bgcolor: "rgba(255,255,255,0.7)",
  bordercolor: "#d1d5db",
  borderwidth: 1
};
const SPECIES_PALETTE = [
  "#2563eb",
  "#dc2626",
  "#16a34a",
  "#9333ea",
  "#ea580c",
  "#0891b2",
  "#be123c",
  "#4d7c0f",
  "#7c3aed",
  "#0f766e"
];

function renderPiCurve(payload) {
  const curve = payload.pi_curve || {positions: [], pi: []};
  const piByPosition = new Map(curve.positions.map((position, index) => [position, curve.pi[index] || 0]));
  const snpY = (payload.snp_positions || []).map((position) => piByPosition.get(position) || 0);
  const indelY = (payload.indel_positions || []).map((position) => piByPosition.get(position) || 0);
  const traces = [
    {
      x: curve.positions,
      y: curve.pi,
      type: "scatter",
      mode: "lines",
      name: "Pi",
      line: {width: 2}
    },
    {
      x: payload.snp_positions || [],
      y: snpY,
      type: "scatter",
      mode: "markers",
      name: "SNP",
      marker: {color: "#dc2626", size: 7}
    },
    {
      x: payload.indel_positions || [],
      y: indelY,
      type: "scatter",
      mode: "markers",
      name: "Indel",
      marker: {color: "#2563eb", size: 8, symbol: "triangle-up"}
    }
  ];
  Plotly.newPlot("pi-curve", traces, {
    margin: PLOT_MARGIN,
    font: PLOT_FONT,
    xaxis: {title: "Position", gridcolor: GRID_COLOR, showgrid: true, zerolinecolor: GRID_COLOR},
    yaxis: {title: "Pi", rangemode: "tozero", gridcolor: GRID_COLOR, showgrid: true, zerolinecolor: GRID_COLOR},
    legend: LEGEND
  }, {responsive: true, displaylogo: false});
}

function renderHaplotypeNetwork(payload) {
  const network = payload.haplotype_network || {nodes: [], edges: []};
  const nodeById = new Map(network.nodes.map((node) => [node.id, node]));
  const shapes = (network.edges || []).map((edge) => {
    const source = nodeById.get(edge.source);
    const target = nodeById.get(edge.target);
    if (!source || !target) {
      return null;
    }
    return {
      type: "line",
      x0: source.x,
      y0: source.y,
      x1: target.x,
      y1: target.y,
      line: {color: "#94a3b8", width: Math.max(1, Math.min(edge.weight, 6))}
    };
  }).filter(Boolean);
  const groups = bySpecies(network.nodes || [], (node) => (node.species || ["Unknown"]).join(", "));
  const traces = Object.entries(groups).map(([species, nodes], index) => ({
    x: nodes.map((node) => node.x),
    y: nodes.map((node) => node.y),
    text: nodes.map((node) => `${node.id}<br>frequency: ${node.frequency}<br>${(node.species || []).join(", ")}`),
    type: "scatter",
    mode: "markers+text",
    name: species,
    textposition: "top center",
    hoverinfo: "text",
    marker: {
      size: nodes.map((node) => 12 + Math.sqrt(node.frequency) * 10),
      color: SPECIES_PALETTE[index % SPECIES_PALETTE.length],
      line: {color: "#ffffff", width: 1}
    }
  }));
  Plotly.newPlot("haplotype-net", traces, {
    margin: PLOT_MARGIN,
    font: PLOT_FONT,
    shapes,
    xaxis: {visible: false, zeroline: false, gridcolor: GRID_COLOR, showgrid: true},
    yaxis: {visible: false, zeroline: false, gridcolor: GRID_COLOR, showgrid: true},
    legend: LEGEND
  }, {responsive: true, displaylogo: false});
}

function renderSpeciesPca(payload) {
  const pca = payload.species_pca || {samples: [], explained_variance: [0, 0]};
  const groups = bySpecies(pca.samples || [], (sample) => sample.species);
  const traces = Object.entries(groups).map(([species, samples], index) => ({
    x: samples.map((sample) => sample.pc1),
    y: samples.map((sample) => sample.pc2),
    text: samples.map((sample) => `${sample.sample_name}<br>${sample.species || "Unknown"}`),
    type: "scatter",
    mode: "markers",
    name: species,
    hovertemplate: "%{text}<br>PC1=%{x:.4f}<br>PC2=%{y:.4f}<extra></extra>",
    marker: {
      size: 10,
      color: SPECIES_PALETTE[index % SPECIES_PALETTE.length],
      line: {color: "#ffffff", width: 1}
    }
  }));
  const variance = pca.explained_variance || [0, 0];
  Plotly.newPlot("species-pca", traces, {
    margin: PLOT_MARGIN,
    font: PLOT_FONT,
    xaxis: {title: `PC1 (${((variance[0] || 0) * 100).toFixed(1)}%)`, gridcolor: GRID_COLOR, showgrid: true, zerolinecolor: GRID_COLOR},
    yaxis: {title: `PC2 (${((variance[1] || 0) * 100).toFixed(1)}%)`, gridcolor: GRID_COLOR, showgrid: true, zerolinecolor: GRID_COLOR},
    legend: LEGEND
  }, {responsive: true, displaylogo: false});
}

fetch("./data.json")
  .then((response) => {
    if (!response.ok) {
      throw new Error("Feature data could not be loaded.");
    }
    return response.json();
  })
  .then((payload) => {
    renderPiCurve(payload);
    renderHaplotypeNetwork(payload);
    renderSpeciesPca(payload);
  });
