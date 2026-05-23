import { useEffect, useState } from 'react';

const repoOwner = 'PI-Prasaad-Krishna';
const repoName = 'FormulaOne-replayer';
const issuesUrl = `https://github.com/${repoOwner}/${repoName}/issues`;
const feedbackUrl = `https://github.com/${repoOwner}/${repoName}/discussions`;
const releaseBaseUrl = `https://github.com/${repoOwner}/${repoName}/releases/latest`;

const features = [
  {
    title: 'Immersive Real-Time Telemetry',
    text: 'Experience race data live. Track speed, gears, and braking with precise visualization widgets.',
  },
  {
    title: 'Cross-Platform Performance',
    text: 'Engineered natively for Windows, macOS, and Linux to deliver butter-smooth performance and low latency.',
  },
  {
    title: 'Professional Interface',
    text: 'Built with a serious, aerodynamic interface. No gimmicks—just pure, data-driven visualization.',
  },
];

function App() {
  const [mounted, setMounted] = useState(false);
  const [releaseInfo, setReleaseInfo] = useState(null);

  useEffect(() => {
    setMounted(true);

    // Fetch GitHub Release Data
    fetch(`https://api.github.com/repos/${repoOwner}/${repoName}/releases/latest`)
      .then((res) => res.json())
      .then((data) => {
        if (data && data.tag_name) {
          setReleaseInfo(data);
        }
      })
      .catch((err) => console.error('Failed to fetch release data:', err));

    // Scroll Animation Observer
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            entry.target.classList.add('is-visible');
          }
        });
      },
      {
        threshold: 0.2,
        rootMargin: '0px 0px -10% 0px',
      },
    );

    document.querySelectorAll('[data-reveal]').forEach((node) => observer.observe(node));

    return () => observer.disconnect();
  }, []);

  const getDownloadUrl = (extension) => {
    if (!releaseInfo || !releaseInfo.assets) return releaseBaseUrl;
    const asset = releaseInfo.assets.find(a => a.name.endsWith(extension));
    return asset ? asset.browser_download_url : releaseBaseUrl;
  };

  const getDownloadSize = (extension) => {
    if (!releaseInfo || !releaseInfo.assets) return 'Checking...';
    const asset = releaseInfo.assets.find(a => a.name.endsWith(extension));
    return asset ? `${(asset.size / 1024 / 1024).toFixed(1)} MB` : 'Available';
  };

  const latestVersion = releaseInfo ? releaseInfo.tag_name : 'Latest';

  const downloads = [
    {
      name: 'Windows',
      filename: `F1 Visualizer Setup (${latestVersion}).msi`,
      details: `Native setup for Windows 10 & 11. Includes automatic updates and desktop shortcuts. Size: ${getDownloadSize('.msi')}`,
      badge: 'Recommended',
      href: getDownloadUrl('.msi'),
    },
    {
      name: 'macOS',
      filename: `F1 Visualizer (${latestVersion}).dmg`,
      details: `Universal build for Apple Silicon and Intel Macs. Clean installation and optimized performance. Size: ${getDownloadSize('.dmg')}`,
      badge: 'Universal',
      href: getDownloadUrl('.dmg'),
    },
    {
      name: 'Linux',
      filename: `F1 Visualizer (${latestVersion}).AppImage`,
      details: `Portable AppImage for modern Linux distributions. No installation required. Size: ${getDownloadSize('.AppImage')}`,
      badge: 'Portable',
      href: getDownloadUrl('.AppImage'),
    },
  ];

  const highlights = [
    {
      title: 'Current Build',
      value: latestVersion,
    },
    {
      title: 'Supported OS',
      value: 'Win, Mac, Lin',
    },
    {
      title: 'Data Integration',
      value: 'Live Telemetry',
    },
  ];

  return (
    <div className={`site ${mounted ? 'is-ready' : ''}`}>
      <div className="grid-overlay" />

      <header className="topbar">
        <a className="brand" href="#top" aria-label="F1 Visualizer home">
          <span className="brand-mark" aria-hidden="true">
            <svg viewBox="0 0 48 48" role="img" aria-hidden="true">
              <path d="M7 27h18l-3 7h9l10-17h-7l4-8H21l-8 10H7z" />
            </svg>
          </span>
          <span>
            <strong>F1 Visualizer</strong>
            <small>Telemetry Dashboard</small>
          </span>
        </a>

        <nav className="nav">
          <a href="#downloads">Download</a>
          <a href="#features">Specifications</a>
          <a href={issuesUrl} target="_blank" rel="noreferrer" title="Report issues">Issues</a>
          <a href={feedbackUrl} target="_blank" rel="noreferrer" title="Provide feedback">Feedback</a>
        </nav>

        <a className="release-pill" href={releaseBaseUrl} target="_blank" rel="noreferrer">
          {latestVersion} Release
        </a>
      </header>

      <main id="top">
        <section className="hero">
          <div className="hero-copy" data-reveal>
            <div className="eyebrow">Professional Race Telemetry</div>
            <h1>The ultimate F1 data visualization suite.</h1>
            <p>
              Download the cross-platform telemetry dashboard for Windows, macOS, and Linux. Built for deep strategy analysis and real-time head-to-head comparisons.
            </p>

            <div className="cta-row">
              <a className="primary-cta" href="#downloads">
                Select Platform
              </a>
              <a className="secondary-cta" href={releaseBaseUrl} target="_blank" rel="noreferrer">
                View Source Repository
              </a>
            </div>

            <div className="stat-row">
              {highlights.map((item) => (
                <div key={item.title} className="stat-card">
                  <span>{item.title}</span>
                  <strong>{item.value}</strong>
                </div>
              ))}
            </div>
          </div>

          <div className="hero-visual" data-reveal>
            <div className="track-panel">
              <div className="track-label">System Dashboard Overview</div>
              <div className="hero-accent-line" aria-hidden="true" />
              <img src="/hero-bg.png" className="hero-image" alt="Dashboard interface preview" />

              <div className="build-card">
                <div>
                  <span className="muted">Build Architecture</span>
                  <strong>Cross-Platform Native</strong>
                </div>
                <div className="build-meta">
                  <span>{latestVersion}</span>
                  <span>Stable Release</span>
                </div>
              </div>
            </div>
          </div>
        </section>

        <section className="downloads-section" id="downloads">
          <div className="section-heading" data-reveal>
            <span>Installation Packages</span>
            <h2>Select your operating system.</h2>
          </div>

          <div className="download-grid">
            {downloads.map((item) => (
              <a
                key={item.name}
                className="download-card"
                href={item.href}
                target="_blank"
                rel="noreferrer"
                data-reveal
              >
                <div className="download-card-top">
                  <span className="download-badge">{item.badge}</span>
                  <span className="download-platform">{item.name}</span>
                </div>

                <div className="download-main">
                  <strong>{item.filename}</strong>
                  <p>{item.details}</p>
                </div>

                <div className="download-footer">
                  <span>Download Package</span>
                  <span aria-hidden="true">
                    <svg viewBox="0 0 24 24">
                      <path d="M13 5l7 7-7 7M20 12H4" />
                    </svg>
                  </span>
                </div>
              </a>
            ))}
          </div>
        </section>

        <section className="features-panel" id="features" data-reveal>
          <div className="section-heading">
            <span className="eyebrow">Technical Specifications</span>
            <h2>Engineered for precise analysis.</h2>
          </div>

          <div className="step-grid">
            {features.map((feature) => (
              <div key={feature.title} className="step-card">
                <strong>{feature.title}</strong>
                <p>{feature.text}</p>
              </div>
            ))}
          </div>
        </section>
      </main>
    </div>
  );
}

export default App;